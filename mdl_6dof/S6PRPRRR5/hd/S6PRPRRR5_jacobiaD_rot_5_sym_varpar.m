% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:11
% EndTime: 2019-02-26 19:56:11
% DurationCPUTime: 0.36s
% Computational Cost: add. (1083->57), mult. (2577->135), div. (441->14), fcn. (3313->11), ass. (0->71)
t169 = sin(pkin(11));
t171 = cos(pkin(11));
t174 = cos(qJ(2));
t172 = cos(pkin(6));
t173 = sin(qJ(2));
t190 = t172 * t173;
t154 = t169 * t174 + t171 * t190;
t170 = sin(pkin(6));
t191 = t170 * t173;
t144 = atan2(-t154, t191);
t140 = sin(t144);
t141 = cos(t144);
t127 = -t140 * t154 + t141 * t191;
t124 = 0.1e1 / t127;
t189 = t172 * t174;
t156 = t169 * t189 + t171 * t173;
t168 = qJ(4) + qJ(5);
t160 = sin(t168);
t161 = cos(t168);
t192 = t169 * t170;
t139 = t156 * t160 + t161 * t192;
t135 = 0.1e1 / t139;
t165 = 0.1e1 / t173;
t125 = 0.1e1 / t127 ^ 2;
t136 = 0.1e1 / t139 ^ 2;
t166 = 0.1e1 / t173 ^ 2;
t157 = -t169 * t190 + t171 * t174;
t201 = 0.2e1 * t157;
t185 = t171 * t189;
t188 = qJD(2) * t173;
t147 = -qJD(2) * t185 + t169 * t188;
t193 = t166 * t174;
t186 = t154 * t193;
t151 = t154 ^ 2;
t163 = 0.1e1 / t170 ^ 2;
t145 = t151 * t163 * t166 + 0.1e1;
t142 = 0.1e1 / t145;
t162 = 0.1e1 / t170;
t195 = t142 * t162;
t119 = (qJD(2) * t186 + t147 * t165) * t195;
t182 = -t140 * t191 - t141 * t154;
t187 = t141 * t170 * t174;
t116 = qJD(2) * t187 + t182 * t119 + t140 * t147;
t200 = t116 * t124 * t125;
t199 = t125 * t157;
t138 = -t156 * t161 + t160 * t192;
t134 = t138 ^ 2;
t131 = t134 * t136 + 0.1e1;
t150 = t157 * qJD(2);
t164 = qJD(4) + qJD(5);
t184 = t164 * t192 - t150;
t194 = t156 * t164;
t133 = t160 * t194 + t184 * t161;
t196 = t136 * t138;
t132 = -t184 * t160 + t161 * t194;
t197 = t132 * t135 * t136;
t198 = 0.1e1 / t131 ^ 2 * (t133 * t196 - t134 * t197);
t183 = t135 * t161 + t160 * t196;
t153 = t169 * t173 - t185;
t181 = t153 * t165 + t186;
t167 = t165 * t166;
t152 = t157 ^ 2;
t149 = t156 * qJD(2);
t148 = t154 * qJD(2);
t128 = 0.1e1 / t131;
t123 = t125 * t152 + 0.1e1;
t120 = t181 * t195;
t117 = t182 * t120 + t140 * t153 + t187;
t115 = (-0.2e1 * t181 / t145 ^ 2 * (-qJD(2) * t151 * t167 * t174 - t147 * t154 * t166) * t163 + (-t147 * t193 + t148 * t165 + (-t153 * t193 + (-0.2e1 * t167 * t174 ^ 2 - t165) * t154) * qJD(2)) * t142) * t162;
t113 = -0.2e1 * t198 + 0.2e1 * (t128 * t133 * t136 + (-t128 * t197 - t136 * t198) * t138) * t138;
t1 = [0, t115, 0, 0, 0, 0; 0, 0.2e1 * (t117 * t199 + t124 * t156) / t123 ^ 2 * (-t149 * t199 - t152 * t200) + (t117 * t200 * t201 - t150 * t124 + (t156 * t116 + t117 * t149 + (-(-t170 * t188 - t115 * t154 + t120 * t147 + (-t120 * t191 + t153) * t119) * t141 - (t119 * t120 * t154 + t148 + (-t115 * t173 + (-qJD(2) * t120 - t119) * t174) * t170) * t140) * t157) * t125) / t123, 0, 0, 0, 0; 0, t183 * t198 * t201 + (t183 * t149 + ((t135 * t164 + 0.2e1 * t138 * t197) * t160 + (-t133 * t160 + (-t138 * t164 + t132) * t161) * t136) * t157) * t128, 0, t113, t113, 0;];
JaD_rot  = t1;
