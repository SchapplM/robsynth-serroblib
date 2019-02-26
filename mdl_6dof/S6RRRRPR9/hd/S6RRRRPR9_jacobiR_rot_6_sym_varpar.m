% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:13
% EndTime: 2019-02-26 22:35:13
% DurationCPUTime: 0.15s
% Computational Cost: add. (205->36), mult. (270->62), div. (0->0), fcn. (395->10), ass. (0->44)
t173 = cos(pkin(6));
t174 = sin(qJ(2));
t177 = cos(qJ(1));
t179 = t177 * t174;
t175 = sin(qJ(1));
t176 = cos(qJ(2));
t180 = t175 * t176;
t161 = t173 * t179 + t180;
t171 = qJ(3) + qJ(4);
t168 = sin(t171);
t169 = cos(t171);
t172 = sin(pkin(6));
t182 = t172 * t177;
t154 = -t161 * t169 + t168 * t182;
t178 = t177 * t176;
t181 = t175 * t174;
t160 = -t173 * t178 + t181;
t170 = pkin(12) + qJ(6);
t166 = sin(t170);
t167 = cos(t170);
t195 = t154 * t166 + t160 * t167;
t194 = t154 * t167 - t160 * t166;
t152 = -t161 * t168 - t169 * t182;
t193 = t152 * t166;
t163 = -t173 * t181 + t178;
t184 = t172 * t175;
t155 = t163 * t168 - t169 * t184;
t192 = t155 * t166;
t185 = t172 * t174;
t158 = -t168 * t185 + t173 * t169;
t191 = t158 * t166;
t188 = t166 * t169;
t187 = t167 * t169;
t186 = t169 * t176;
t183 = t172 * t176;
t162 = t173 * t180 + t179;
t159 = t173 * t168 + t169 * t185;
t157 = t158 * t167;
t156 = t163 * t169 + t168 * t184;
t151 = t155 * t167;
t150 = t152 * t167;
t149 = t156 * t167 + t162 * t166;
t148 = -t156 * t166 + t162 * t167;
t1 = [t194, -t162 * t187 + t163 * t166, -t151, -t151, 0, t148; t149, -t160 * t187 + t161 * t166, t150, t150, 0, t195; 0 (t166 * t174 + t167 * t186) * t172, t157, t157, 0, -t159 * t166 - t167 * t183; -t195, t162 * t188 + t163 * t167, t192, t192, 0, -t149; t148, t160 * t188 + t161 * t167, -t193, -t193, 0, t194; 0 (-t166 * t186 + t167 * t174) * t172, -t191, -t191, 0, -t159 * t167 + t166 * t183; t152, -t162 * t168, t156, t156, 0, 0; t155, -t160 * t168, -t154, -t154, 0, 0; 0, t168 * t183, t159, t159, 0, 0;];
JR_rot  = t1;
