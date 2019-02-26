% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRPRR13_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:43
% EndTime: 2019-02-26 20:55:44
% DurationCPUTime: 0.05s
% Computational Cost: add. (33->20), mult. (98->45), div. (0->0), fcn. (140->12), ass. (0->30)
t160 = sin(pkin(7));
t164 = cos(pkin(6));
t180 = t160 * t164;
t161 = sin(pkin(6));
t167 = sin(qJ(1));
t179 = t161 * t167;
t170 = cos(qJ(1));
t178 = t161 * t170;
t162 = cos(pkin(12));
t163 = cos(pkin(7));
t177 = t162 * t163;
t159 = sin(pkin(12));
t176 = t167 * t159;
t175 = t167 * t162;
t174 = t170 * t159;
t173 = t170 * t162;
t155 = t164 * t173 - t176;
t172 = -t155 * t163 + t160 * t178;
t157 = -t164 * t175 - t174;
t171 = t157 * t163 + t160 * t179;
t169 = cos(qJ(3));
t168 = cos(qJ(5));
t166 = sin(qJ(3));
t165 = sin(qJ(5));
t158 = -t164 * t176 + t173;
t156 = t164 * t174 + t175;
t154 = -t161 * t162 * t160 + t164 * t163;
t153 = -t157 * t160 + t163 * t179;
t152 = -t155 * t160 - t163 * t178;
t1 = [0, 0, t153, 0, t158 * t169 + t171 * t166, t153 * t165 - (t158 * t166 - t171 * t169) * t168; 0, 0, t152, 0, t156 * t169 - t172 * t166, t152 * t165 - (t156 * t166 + t172 * t169) * t168; 1, 0, t154, 0, t166 * t180 + (t159 * t169 + t166 * t177) * t161, t154 * t165 - (-t169 * t180 + (t159 * t166 - t169 * t177) * t161) * t168;];
Jg_rot  = t1;
