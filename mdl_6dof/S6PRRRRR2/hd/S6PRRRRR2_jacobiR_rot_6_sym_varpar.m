% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:15
% EndTime: 2019-02-26 20:19:15
% DurationCPUTime: 0.11s
% Computational Cost: add. (209->32), mult. (270->64), div. (0->0), fcn. (395->10), ass. (0->44)
t173 = sin(pkin(12));
t175 = cos(pkin(12));
t178 = cos(qJ(2));
t176 = cos(pkin(6));
t177 = sin(qJ(2));
t180 = t176 * t177;
t163 = t173 * t178 + t175 * t180;
t172 = qJ(3) + qJ(4);
t168 = sin(t172);
t170 = cos(t172);
t174 = sin(pkin(6));
t183 = t174 * t175;
t155 = -t163 * t168 - t170 * t183;
t171 = qJ(5) + qJ(6);
t167 = sin(t171);
t190 = t155 * t167;
t165 = -t173 * t180 + t175 * t178;
t184 = t173 * t174;
t157 = -t165 * t168 + t170 * t184;
t189 = t157 * t167;
t182 = t174 * t177;
t160 = -t168 * t182 + t170 * t176;
t188 = t160 * t167;
t187 = t167 * t170;
t169 = cos(t171);
t186 = t169 * t170;
t185 = t170 * t178;
t181 = t174 * t178;
t179 = t176 * t178;
t164 = t173 * t179 + t175 * t177;
t162 = t173 * t177 - t175 * t179;
t161 = t168 * t176 + t170 * t182;
t159 = t160 * t169;
t158 = t165 * t170 + t168 * t184;
t156 = t163 * t170 - t168 * t183;
t154 = -t161 * t169 + t167 * t181;
t153 = -t161 * t167 - t169 * t181;
t152 = t157 * t169;
t151 = t155 * t169;
t150 = -t158 * t169 - t164 * t167;
t149 = -t158 * t167 + t164 * t169;
t148 = -t156 * t169 - t162 * t167;
t147 = -t156 * t167 + t162 * t169;
t1 = [0, -t164 * t186 + t165 * t167, t152, t152, t149, t149; 0, -t162 * t186 + t163 * t167, t151, t151, t147, t147; 0 (t167 * t177 + t169 * t185) * t174, t159, t159, t153, t153; 0, t164 * t187 + t165 * t169, -t189, -t189, t150, t150; 0, t162 * t187 + t163 * t169, -t190, -t190, t148, t148; 0 (-t167 * t185 + t169 * t177) * t174, -t188, -t188, t154, t154; 0, -t164 * t168, t158, t158, 0, 0; 0, -t162 * t168, t156, t156, 0, 0; 0, t168 * t181, t161, t161, 0, 0;];
JR_rot  = t1;
