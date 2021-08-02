#!/usr/bin/env python3

import cgi
import os
import sys

from axle.add import add
from axle.init import init
from axle.pull import pull
from openpyxl import load_workbook
from urllib.parse import parse_qsl


def build_form(args):
    output = [
        "form",
        {"action": "load_sheet.py", "method": "POST", "enctype": "multipart/form-data"},
        # ["p", str(args)],
        ["input", {"type": "hidden", "name": "action", "value": "upload"}],
        build_input(args, "Upload MRO", input_type="file"),
        build_input(args, "Submit", input_type="submit"),
    ]
    return render_output(output)


def build_input(args, label, input_type="text"):
    output = ["div", {"class": "form-group row"}]
    name = label.lower().replace(" ", "_")
    value = args[name] if name in args else ""
    left = "col-md-3"
    right = "col-md-9"
    control = [right, "form-control"]
    if "valid" in args and name in args["valid"]:
        control.append("is-valid")
    elif "invalid" in args and name in args["invalid"]:
        control.append("is-invalid")
    label_classes = " ".join([left, "col-form-label"])
    control_classes = " ".join(control)

    if input_type == "textarea":
        output.append(["label", {"for": name, "class": label_classes}, label])
        output.append(["textarea", {"class": control_classes, "name": name}, value])
    elif input_type == "file":
        output.append(["label", {"for": name, "class": label_classes}, label])
        control_classes = control_classes.replace("form-control", "form-control-file")
        output.append(
            ["input", {"type": "file", "class": control_classes, "name": name, "value": label}]
        )
    elif input_type == "submit":
        output.append(
            ["input", {"type": "submit", "class": "btn btn-primary", "name": name, "value": label}]
        )
    else:
        output.append(["label", {"for": name, "class": label_classes}, label])
        output.append(
            ["input", {"type": "text", "class": control_classes, "name": name, "value": value}]
        )

    if "valid" in args and name in args["valid"] and isinstance(args["valid"][name], str):
        output.append(["div", {"class": left}])
        output.append(["div", {"class": right + " valid-feedback"}, args["valid"][name]])
    if "invalid" in args and name in args["invalid"] and isinstance(args["invalid"][name], str):
        output.append(["div", {"class": left}])
        output.append(["div", {"class": right + " invalid-feedback"}, args["invalid"][name]])

    return output


def render_html(element, depth=0):
    """Render hiccup-style HTML vector as HTML."""
    indent = "  " * depth
    if not isinstance(element, list):
        raise Exception(f"Element is not a list: {element}")
    if len(element) == 0:
        raise Exception(f"Element is an empty list")
    tag = element.pop(0)
    if not isinstance(tag, str):
        raise Exception(f"Tag '{tag}' is not a string in '{element}'")
    output = f"{indent}<{tag}"

    if len(element) > 0 and isinstance(element[0], dict):
        attrs = element.pop(0)
        if tag == "a" and "href" not in attrs and "resource" in attrs:
            attrs["href"] = curie2href(attrs["resource"])
        for key, value in attrs.items():
            if key in ["checked"]:
                if value:
                    output += f" {key}"
            else:
                output += f' {key}="{value}"'

    if tag in ["meta", "link"]:
        output += "/>"
        return output
    output += ">"
    spacing = ""
    if len(element) > 0:
        for child in element:
            if isinstance(child, str):
                output += child
            elif isinstance(child, list):
                try:
                    output += "\n" + render_html(child, depth=depth + 1)
                    spacing = f"\n{indent}"
                except Exception as e:
                    raise Exception(f"Bad child in '{element}'", e)
            else:
                raise Exception(f"Bad type for child '{child}' in '{element}'")
    output += f"{spacing}</{tag}>"
    return output


def render_output(output):
    # Render output
    if os.environ.get("GATEWAY_INTERFACE") == "CGI/1.1":
        print("Content-Type: text/html")
        print("")
        print("<head>")
        print('<link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css">')
        print("</head>")
        print("<body>")
        print('<div class="container" style="padding-top:30px;">')
        print(render_html(output))
        print("</div>")
        print("</body>")
    else:
        print(render_text(output))


def main():
    if os.environ["REQUEST_METHOD"] == "POST":
        fields = cgi.FieldStorage()
        args = {}
        for key in fields.keys():
            args[key] = fields[key].value
    elif "QUERY_STRING" in os.environ:
        args = dict(parse_qsl(os.environ["QUERY_STRING"]))
    else:
        return render_output(["h1", "Error: missing QUERY_STRING or POST request"])

    action = args.get("action")
    if not action:
        return render_output(["h1", "Error: please specify an action"])

    if action == "create":
        # Show form to upload the new sheet
        return build_form(args)
    elif action == "upload":
        # Save XLSX to mro.xlsx
        wb = load_workbook(fields["upload_mro"].file)
        wb.save("../../../mro.xlsx")

        # Create new AXLE project and update the existing TSVs
        os.chdir("../../..")
        if not os.path.exists(".axle"):
            init("mro", filepath="mro.xlsx")
            # Add the sheets in the correct order
            add("index.tsv")
            add("iedb/iedb.tsv")
            for sheet in ["genetic-locus", "haplotype", "serotype", "chain", "chain-sequence", "molecule", "haplotype-molecule", "serotype-molecule", "mutant-molecule", "core", "external", "iedb-manual", "evidence", "rejected"]:
                add(f"ontology/{sheet}.tsv")
        pull()


if __name__ == '__main__':
    main()

