import logging
from Aries.outputs import Print
from .utils import api_collection, api_response
from . import basespace, bs_sample
logger = logging.getLogger(__name__)


def basespace_command(options):
    url = options['url']

    if url:
        logger.debug(api_response(url))
        return

    if not options['collection']:
        print("Argument: collection is required when --url is not specified.")
        return

    collection = options['collection']
    basespace_id = options['basespace_id']

    if basespace_id:
        if collection == 'samples':
            href = 'v1pre3/%s/%s/' % (collection, basespace_id)
            print("Getting data from %s..." % href)
            Print.print(basespace.get_details("samples", basespace_id))
            Print.print(bs_sample.get_files(basespace_id, extension=None))
        else:
            properties = options['properties']
            href = 'v1pre3/%s/%s/properties' % (collection, basespace_id)
            print("Getting data from %s..." % href)
            property_items = api_collection(href)
            available_properties = dict()
            for item in property_items:
                available_properties[str(item.get("Href")).split("/")[-1]] = item

            if not properties:
                details = basespace.get_details(collection, basespace_id)
                details.pop("Properties", None)
                Print.print(details)
                Print.green("Available Properties: ")
                for p, v in available_properties.items():
                    items_count = v.get("ItemsTotalCount")
                    print("%-35s| Items Count: %s" % (p, items_count))
            else:
                for p in properties:
                    Print.green(p)
                    if p in available_properties.keys():
                        Print.print(available_properties[p])
                    else:
                        print("Property not found.")

    else:
        items = basespace.get_list(collection)
        for item in items:
            if collection == "runs":
                print("ID: %-15s | ExperimentName: %s" % (item.get("Id"), item.get("ExperimentName")))
            else:
                print("ID: %-15s | Name: %s" % (item.get("Id"), item.get("Name")))
