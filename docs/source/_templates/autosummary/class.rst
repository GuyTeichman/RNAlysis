{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
   {% block attributes %}
   {% if attributes %}
    .. automethod:: __init__
    {% endif %}
    {% endblock %}

   {% if methods %}
.. autosummary::
    :toctree: .
    {% for item in all_methods %}
    {%- if not item.startswith('_') or item in ['__call__'] %}
    {{ name }}.{{ item }}
    {%- endif -%}
    {%- endfor %}
    {% endif %}
    {% endblock %}
