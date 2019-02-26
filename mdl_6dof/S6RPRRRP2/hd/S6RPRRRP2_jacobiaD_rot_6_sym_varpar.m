% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:08:37
% EndTime: 2019-02-26 21:08:38
% DurationCPUTime: 0.72s
% Computational Cost: add. (2635->97), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
t159 = sin(qJ(3));
t154 = t159 ^ 2;
t160 = cos(qJ(3));
t156 = 0.1e1 / t160 ^ 2;
t199 = t154 * t156;
t152 = qJ(1) + pkin(10);
t147 = sin(t152);
t223 = 0.2e1 * t147;
t222 = t159 * t199;
t148 = cos(t152);
t158 = qJ(4) + qJ(5);
t149 = sin(t158);
t150 = cos(t158);
t201 = t150 * t160;
t135 = t147 * t149 + t148 * t201;
t151 = qJD(4) + qJD(5);
t177 = t151 * t160 - qJD(1);
t194 = qJD(3) * t159;
t221 = t177 * t149 + t150 * t194;
t204 = t147 * t159;
t139 = atan2(-t204, -t160);
t138 = cos(t139);
t137 = sin(t139);
t187 = t137 * t204;
t125 = -t138 * t160 - t187;
t122 = 0.1e1 / t125;
t129 = 0.1e1 / t135;
t155 = 0.1e1 / t160;
t123 = 0.1e1 / t125 ^ 2;
t130 = 0.1e1 / t135 ^ 2;
t220 = -0.2e1 * t159;
t145 = t147 ^ 2;
t143 = t145 * t199 + 0.1e1;
t141 = 0.1e1 / t143;
t219 = t141 - 0.1e1;
t202 = t149 * t160;
t170 = t147 * t202 + t148 * t150;
t183 = t149 * t194;
t113 = t170 * qJD(1) - t135 * t151 + t148 * t183;
t134 = -t147 * t150 + t148 * t202;
t128 = t134 ^ 2;
t121 = t128 * t130 + 0.1e1;
t209 = t130 * t134;
t176 = -qJD(1) * t160 + t151;
t172 = t176 * t150;
t114 = t147 * t172 - t148 * t221;
t216 = t114 * t129 * t130;
t218 = (-t113 * t209 - t128 * t216) / t121 ^ 2;
t196 = qJD(1) * t159;
t184 = t148 * t196;
t193 = qJD(3) * t160;
t195 = qJD(3) * t147;
t115 = (-(-t147 * t193 - t184) * t155 + t195 * t199) * t141;
t207 = t138 * t159;
t109 = (-t115 * t147 + qJD(3)) * t207 + (-t184 + (t115 - t195) * t160) * t137;
t217 = t109 * t122 * t123;
t215 = t115 * t137;
t214 = t115 * t159;
t213 = t123 * t148;
t212 = t123 * t159;
t169 = qJD(3) * (t159 + t222) * t155;
t197 = qJD(1) * t148;
t174 = t147 * t154 * t197;
t211 = (t145 * t169 + t156 * t174) / t143 ^ 2;
t181 = 0.1e1 + t199;
t127 = t181 * t147 * t141;
t210 = t127 * t147;
t208 = t137 * t160;
t146 = t148 ^ 2;
t206 = t146 * t154;
t203 = t148 * t149;
t200 = t154 * t155;
t198 = qJD(1) * t147;
t118 = t123 * t206 + 0.1e1;
t192 = 0.2e1 * (-t206 * t217 + (t146 * t159 * t193 - t174) * t123) / t118 ^ 2;
t191 = 0.2e1 * t218;
t190 = 0.2e1 * t217;
t189 = t134 * t216;
t188 = t148 * t212;
t186 = t141 * t200;
t180 = t159 * t192;
t179 = t211 * t223;
t178 = t211 * t220;
t175 = t147 * t186;
t173 = t181 * t148;
t171 = -t129 * t149 + t150 * t209;
t168 = t171 * t159;
t133 = -t147 * t201 + t203;
t119 = 0.1e1 / t121;
t116 = 0.1e1 / t118;
t112 = (t219 * t159 * t137 - t138 * t175) * t148;
t111 = -t147 * t208 + t207 + (-t138 * t204 + t208) * t127;
t110 = -t181 * t179 + (qJD(1) * t173 + t169 * t223) * t141;
t106 = -0.2e1 * t218 + 0.2e1 * (-t113 * t119 * t130 + (-t119 * t216 - t130 * t218) * t134) * t134;
t1 = [t148 * t155 * t178 + (-t147 * t155 * t196 + qJD(3) * t173) * t141, 0, t110, 0, 0, 0; (t122 * t180 + (-t122 * t193 + (qJD(1) * t112 + t109) * t212) * t116) * t147 + (t123 * t180 * t112 + (-((t115 * t175 + t219 * t193 + t178) * t137 + (t179 * t200 - t214 + (t214 + (t220 - t222) * t195) * t141) * t138) * t188 + (-t123 * t193 + t159 * t190) * t112 + (-t122 + ((-t145 + t146) * t138 * t186 + t219 * t187) * t123) * t196) * t116) * t148, 0 (t111 * t212 - t122 * t160) * t148 * t192 + ((-t122 * t198 + (-qJD(3) * t111 - t109) * t213) * t160 + (-t148 * qJD(3) * t122 - (-t110 * t138 * t147 + t137 * t195 + t210 * t215 - t215 + (-qJD(3) * t137 - t138 * t197) * t127) * t188 + (t123 * t198 + t148 * t190) * t111 - ((t110 - t197) * t137 + ((0.1e1 - t210) * qJD(3) + (t127 - t147) * t115) * t138) * t160 * t213) * t159) * t116, 0, 0, 0; (t129 * t170 + t133 * t209) * t191 + (0.2e1 * t133 * t189 + (t133 * t113 + t170 * t114 + (-t147 * t221 - t148 * t172) * t134) * t130 + (t176 * t203 + (-t177 * t150 + t183) * t147) * t129) * t119, 0, -t148 * t168 * t191 + (-t168 * t198 + (t171 * t193 + ((-t129 * t151 - 0.2e1 * t189) * t150 + (-t113 * t150 + (-t134 * t151 + t114) * t149) * t130) * t159) * t148) * t119, t106, t106, 0;];
JaD_rot  = t1;
