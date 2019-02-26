% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:05
% EndTime: 2019-02-26 20:07:05
% DurationCPUTime: 0.15s
% Computational Cost: add. (130->39), mult. (304->88), div. (0->0), fcn. (425->12), ass. (0->46)
t160 = pkin(13) + qJ(5);
t158 = sin(t160);
t162 = sin(pkin(7));
t186 = t158 * t162;
t159 = cos(t160);
t185 = t159 * t162;
t163 = sin(pkin(6));
t184 = t162 * t163;
t166 = cos(pkin(6));
t183 = t162 * t166;
t165 = cos(pkin(7));
t182 = t163 * t165;
t167 = sin(qJ(3));
t181 = t165 * t167;
t169 = cos(qJ(3));
t180 = t165 * t169;
t168 = sin(qJ(2));
t179 = t166 * t168;
t170 = cos(qJ(2));
t178 = t166 * t170;
t177 = t167 * t168;
t176 = t167 * t170;
t175 = t168 * t169;
t174 = t169 * t170;
t173 = t168 * t184;
t161 = sin(pkin(12));
t164 = cos(pkin(12));
t153 = -t161 * t168 + t164 * t178;
t172 = t153 * t165 - t164 * t184;
t155 = -t161 * t178 - t164 * t168;
t171 = t155 * t165 + t161 * t184;
t156 = -t161 * t179 + t164 * t170;
t154 = t161 * t170 + t164 * t179;
t152 = t166 * t165 - t170 * t184;
t151 = (-t165 * t177 + t174) * t163;
t150 = -t155 * t162 + t161 * t182;
t149 = -t153 * t162 - t164 * t182;
t148 = t167 * t183 + (t165 * t176 + t175) * t163;
t147 = t169 * t183 + (t165 * t174 - t177) * t163;
t146 = t155 * t169 - t156 * t181;
t145 = t153 * t169 - t154 * t181;
t144 = t156 * t169 + t171 * t167;
t143 = -t156 * t167 + t171 * t169;
t142 = t154 * t169 + t172 * t167;
t141 = -t154 * t167 + t172 * t169;
t1 = [0, t146 * t159 + t156 * t186, t143 * t159, 0, -t144 * t158 + t150 * t159, 0; 0, t145 * t159 + t154 * t186, t141 * t159, 0, -t142 * t158 + t149 * t159, 0; 0, t151 * t159 + t158 * t173, t147 * t159, 0, -t148 * t158 + t152 * t159, 0; 0, -t146 * t158 + t156 * t185, -t143 * t158, 0, -t144 * t159 - t150 * t158, 0; 0, -t145 * t158 + t154 * t185, -t141 * t158, 0, -t142 * t159 - t149 * t158, 0; 0, -t151 * t158 + t159 * t173, -t147 * t158, 0, -t148 * t159 - t152 * t158, 0; 0, t155 * t167 + t156 * t180, t144, 0, 0, 0; 0, t153 * t167 + t154 * t180, t142, 0, 0, 0; 0 (t165 * t175 + t176) * t163, t148, 0, 0, 0;];
JR_rot  = t1;
