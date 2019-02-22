% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:41
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR7_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiR_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:41:41
% EndTime: 2019-02-22 09:41:41
% DurationCPUTime: 0.12s
% Computational Cost: add. (118->44), mult. (365->102), div. (0->0), fcn. (498->14), ass. (0->49)
t152 = sin(pkin(14));
t160 = cos(pkin(7));
t186 = t152 * t160;
t154 = sin(pkin(8));
t155 = sin(pkin(7));
t185 = t154 * t155;
t156 = sin(pkin(6));
t184 = t155 * t156;
t159 = cos(pkin(8));
t183 = t155 * t159;
t161 = cos(pkin(6));
t182 = t155 * t161;
t181 = t156 * t160;
t163 = sin(qJ(2));
t180 = t156 * t163;
t157 = cos(pkin(14));
t179 = t157 * t160;
t178 = t160 * t163;
t165 = cos(qJ(2));
t177 = t160 * t165;
t176 = t161 * t163;
t175 = t161 * t165;
t174 = t155 * t180;
t153 = sin(pkin(13));
t158 = cos(pkin(13));
t147 = -t153 * t163 + t158 * t175;
t148 = t153 * t165 + t158 * t176;
t168 = t147 * t160 - t158 * t184;
t173 = (-t148 * t152 + t168 * t157) * t159 + (-t147 * t155 - t158 * t181) * t154;
t149 = -t153 * t175 - t158 * t163;
t150 = -t153 * t176 + t158 * t165;
t167 = t149 * t160 + t153 * t184;
t172 = (-t150 * t152 + t167 * t157) * t159 + (-t149 * t155 + t153 * t181) * t154;
t171 = (t157 * t182 + (-t152 * t163 + t157 * t177) * t156) * t159 + (t161 * t160 - t165 * t184) * t154;
t136 = -t147 * t152 - t148 * t179;
t170 = t136 * t159 + t148 * t185;
t138 = -t149 * t152 - t150 * t179;
t169 = t138 * t159 + t150 * t185;
t144 = (-t152 * t165 - t157 * t178) * t156;
t166 = t144 * t159 + t154 * t174;
t164 = cos(qJ(4));
t162 = sin(qJ(4));
t145 = (-t152 * t178 + t157 * t165) * t156;
t141 = t157 * t180 + (t156 * t177 + t182) * t152;
t139 = t149 * t157 - t150 * t186;
t137 = t147 * t157 - t148 * t186;
t135 = t150 * t157 + t167 * t152;
t133 = t148 * t157 + t168 * t152;
t1 = [0, t139 * t164 + t169 * t162, 0, -t135 * t162 + t172 * t164, 0, 0; 0, t137 * t164 + t170 * t162, 0, -t133 * t162 + t173 * t164, 0, 0; 0, t145 * t164 + t166 * t162, 0, -t141 * t162 + t171 * t164, 0, 0; 0, -t139 * t162 + t169 * t164, 0, -t135 * t164 - t172 * t162, 0, 0; 0, -t137 * t162 + t170 * t164, 0, -t133 * t164 - t173 * t162, 0, 0; 0, -t145 * t162 + t166 * t164, 0, -t141 * t164 - t171 * t162, 0, 0; 0, -t138 * t154 + t150 * t183, 0, 0, 0, 0; 0, -t136 * t154 + t148 * t183, 0, 0, 0, 0; 0, -t144 * t154 + t159 * t174, 0, 0, 0, 0;];
JR_rot  = t1;
