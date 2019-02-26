% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP8_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:49
% EndTime: 2019-02-26 22:43:49
% DurationCPUTime: 0.16s
% Computational Cost: add. (161->35), mult. (270->63), div. (0->0), fcn. (395->10), ass. (0->43)
t165 = cos(pkin(6));
t167 = sin(qJ(2));
t171 = cos(qJ(1));
t173 = t171 * t167;
t168 = sin(qJ(1));
t170 = cos(qJ(2));
t175 = t168 * t170;
t156 = t165 * t173 + t175;
t163 = qJ(3) + qJ(4);
t161 = sin(t163);
t162 = cos(t163);
t164 = sin(pkin(6));
t178 = t164 * t171;
t149 = -t156 * t162 + t161 * t178;
t172 = t171 * t170;
t176 = t168 * t167;
t155 = -t165 * t172 + t176;
t166 = sin(qJ(5));
t169 = cos(qJ(5));
t189 = t149 * t166 + t155 * t169;
t188 = t149 * t169 - t155 * t166;
t147 = -t156 * t161 - t162 * t178;
t187 = t147 * t166;
t158 = -t165 * t176 + t172;
t179 = t164 * t168;
t150 = t158 * t161 - t162 * t179;
t186 = t150 * t166;
t180 = t164 * t167;
t153 = -t161 * t180 + t165 * t162;
t185 = t153 * t166;
t182 = t162 * t166;
t181 = t162 * t169;
t177 = t166 * t170;
t174 = t169 * t170;
t157 = t165 * t175 + t173;
t154 = t165 * t161 + t162 * t180;
t152 = t153 * t169;
t151 = t158 * t162 + t161 * t179;
t146 = t150 * t169;
t145 = t147 * t169;
t144 = t151 * t169 + t157 * t166;
t143 = -t151 * t166 + t157 * t169;
t1 = [t188, -t157 * t181 + t158 * t166, -t146, -t146, t143, 0; t144, -t155 * t181 + t156 * t166, t145, t145, t189, 0; 0 (t162 * t174 + t166 * t167) * t164, t152, t152, -t154 * t166 - t164 * t174, 0; -t189, t157 * t182 + t158 * t169, t186, t186, -t144, 0; t143, t155 * t182 + t156 * t169, -t187, -t187, t188, 0; 0 (-t162 * t177 + t167 * t169) * t164, -t185, -t185, -t154 * t169 + t164 * t177, 0; t147, -t157 * t161, t151, t151, 0, 0; t150, -t155 * t161, -t149, -t149, 0, 0; 0, t164 * t170 * t161, t154, t154, 0, 0;];
JR_rot  = t1;
