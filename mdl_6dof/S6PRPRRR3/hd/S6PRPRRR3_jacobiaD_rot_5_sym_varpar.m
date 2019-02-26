% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:48
% EndTime: 2019-02-26 19:54:49
% DurationCPUTime: 0.39s
% Computational Cost: add. (1419->57), mult. (2577->133), div. (441->14), fcn. (3313->11), ass. (0->69)
t173 = sin(qJ(2));
t174 = cos(qJ(2));
t171 = sin(pkin(11));
t204 = cos(pkin(6));
t189 = t171 * t204;
t203 = cos(pkin(11));
t158 = -t173 * t189 + t203 * t174;
t185 = t204 * t203;
t154 = t171 * t173 - t174 * t185;
t172 = sin(pkin(6));
t193 = t172 * t174;
t144 = atan2(-t154, -t193);
t142 = sin(t144);
t143 = cos(t144);
t129 = -t142 * t154 - t143 * t193;
t126 = 0.1e1 / t129;
t164 = pkin(12) + qJ(4) + qJ(5);
t162 = sin(t164);
t163 = cos(t164);
t194 = t171 * t172;
t141 = t158 * t163 + t162 * t194;
t137 = 0.1e1 / t141;
t168 = 0.1e1 / t174;
t127 = 0.1e1 / t129 ^ 2;
t138 = 0.1e1 / t141 ^ 2;
t169 = 0.1e1 / t174 ^ 2;
t156 = t171 * t174 + t173 * t185;
t149 = t156 * qJD(2);
t195 = t169 * t173;
t190 = t154 * t195;
t152 = t154 ^ 2;
t166 = 0.1e1 / t172 ^ 2;
t147 = t152 * t166 * t169 + 0.1e1;
t145 = 0.1e1 / t147;
t165 = 0.1e1 / t172;
t197 = t145 * t165;
t121 = (qJD(2) * t190 + t149 * t168) * t197;
t183 = t142 * t193 - t143 * t154;
t191 = t143 * t172 * t173;
t118 = qJD(2) * t191 + t121 * t183 - t142 * t149;
t202 = t118 * t126 * t127;
t140 = t158 * t162 - t163 * t194;
t136 = t140 ^ 2;
t132 = t136 * t138 + 0.1e1;
t181 = -t173 * t203 - t174 * t189;
t150 = t181 * qJD(2);
t167 = qJD(4) + qJD(5);
t186 = t167 * t194 + t150;
t196 = t158 * t167;
t133 = t162 * t186 + t163 * t196;
t198 = t138 * t140;
t134 = -t162 * t196 + t163 * t186;
t199 = t134 * t137 * t138;
t201 = (t133 * t198 - t136 * t199) / t132 ^ 2;
t200 = t127 * t181;
t192 = -0.2e1 * t201;
t184 = -t137 * t162 + t163 * t198;
t182 = t156 * t168 + t190;
t170 = t168 * t169;
t153 = t181 ^ 2;
t151 = t158 * qJD(2);
t148 = qJD(2) * t154;
t130 = 0.1e1 / t132;
t125 = t153 * t127 + 0.1e1;
t122 = t182 * t197;
t119 = t122 * t183 - t142 * t156 + t191;
t117 = (-0.2e1 * t182 / t147 ^ 2 * (qJD(2) * t152 * t170 * t173 + t149 * t154 * t169) * t166 + (t149 * t195 - t148 * t168 + (t156 * t195 + (0.2e1 * t170 * t173 ^ 2 + t168) * t154) * qJD(2)) * t145) * t165;
t115 = t192 + 0.2e1 * (t130 * t133 * t138 + (-t130 * t199 - t138 * t201) * t140) * t140;
t1 = [0, t117, 0, 0, 0, 0; 0, 0.2e1 * (-t119 * t200 - t126 * t158) / t125 ^ 2 * (-t151 * t200 - t153 * t202) + (t150 * t126 + (-t158 * t118 - t119 * t151) * t127 - (0.2e1 * t119 * t202 + (-(qJD(2) * t193 - t117 * t154 - t122 * t149 + (t122 * t193 - t156) * t121) * t143 - (t121 * t122 * t154 + t148 + (t117 * t174 + (-qJD(2) * t122 - t121) * t173) * t172) * t142) * t127) * t181) / t125, 0, 0, 0, 0; 0, -t184 * t181 * t192 + (t184 * t151 - ((-t137 * t167 - 0.2e1 * t140 * t199) * t163 + (t133 * t163 + (-t140 * t167 + t134) * t162) * t138) * t181) * t130, 0, t115, t115, 0;];
JaD_rot  = t1;
