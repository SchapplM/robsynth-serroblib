% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:37:36
% EndTime: 2019-02-26 20:37:36
% DurationCPUTime: 0.72s
% Computational Cost: add. (1192->94), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->97)
t154 = sin(qJ(4));
t147 = 0.1e1 / t154 ^ 2;
t156 = cos(qJ(4));
t151 = t156 ^ 2;
t201 = t147 * t151;
t219 = t156 * t201;
t155 = sin(qJ(1));
t177 = 0.1e1 + t201;
t218 = t155 * t177;
t149 = t155 ^ 2;
t141 = t149 * t201 + 0.1e1;
t138 = 0.1e1 / t141;
t146 = 0.1e1 / t154;
t157 = cos(qJ(1));
t191 = qJD(1) * t157;
t179 = t156 * t191;
t189 = qJD(4) * t155;
t113 = ((-t154 * t189 + t179) * t146 - t189 * t201) * t138;
t217 = -t113 - t189;
t195 = t155 * t156;
t140 = atan2(t195, t154);
t137 = cos(t140);
t136 = sin(t140);
t182 = t136 * t195;
t123 = t137 * t154 + t182;
t120 = 0.1e1 / t123;
t153 = qJ(5) + qJ(6);
t144 = cos(t153);
t198 = t154 * t157;
t180 = t144 * t198;
t143 = sin(t153);
t197 = t155 * t143;
t133 = t180 - t197;
t127 = 0.1e1 / t133;
t121 = 0.1e1 / t123 ^ 2;
t128 = 0.1e1 / t133 ^ 2;
t216 = -0.2e1 * t156;
t215 = t138 - 0.1e1;
t152 = t157 ^ 2;
t200 = t151 * t152;
t118 = t121 * t200 + 0.1e1;
t199 = t151 * t155;
t170 = t191 * t199;
t188 = qJD(4) * t156;
t203 = t137 * t156;
t107 = (t113 * t155 + qJD(4)) * t203 + (t154 * t217 + t179) * t136;
t212 = t107 * t120 * t121;
t214 = (-t200 * t212 + (-t152 * t154 * t188 - t170) * t121) / t118 ^ 2;
t145 = qJD(5) + qJD(6);
t172 = qJD(1) * t154 + t145;
t187 = qJD(4) * t157;
t178 = t156 * t187;
t111 = -t143 * t178 - t144 * t191 - t145 * t180 + t172 * t197;
t196 = t155 * t144;
t132 = t143 * t198 + t196;
t126 = t132 ^ 2;
t119 = t126 * t128 + 0.1e1;
t206 = t128 * t132;
t173 = t145 * t154 + qJD(1);
t202 = t143 * t157;
t112 = -t173 * t202 + (-t172 * t155 + t178) * t144;
t211 = t112 * t127 * t128;
t213 = (-t111 * t206 - t126 * t211) / t119 ^ 2;
t210 = t113 * t156;
t209 = t121 * t156;
t168 = (t156 + t219) * t146;
t208 = (-t168 * t149 * qJD(4) + t147 * t170) / t141 ^ 2;
t207 = t127 * t143;
t205 = t132 * t144;
t204 = t136 * t155;
t194 = t156 * t157;
t193 = qJD(1) * t155;
t192 = qJD(1) * t156;
t190 = qJD(4) * t154;
t186 = 0.2e1 * t213;
t185 = -0.2e1 * t212;
t184 = 0.2e1 * t208;
t183 = t121 * t194;
t181 = t138 * t146 * t151;
t176 = t214 * t216;
t175 = 0.2e1 * t132 * t211;
t174 = -0.2e1 * t146 * t208;
t171 = t155 * t181;
t169 = t177 * t157;
t167 = t128 * t205 - t207;
t166 = t167 * t157;
t165 = -t155 * t188 - t172 * t157;
t131 = -t154 * t196 - t202;
t130 = t144 * t157 - t154 * t197;
t125 = t138 * t218;
t116 = 0.1e1 / t119;
t114 = 0.1e1 / t118;
t110 = (-t215 * t156 * t136 + t137 * t171) * t157;
t109 = -t154 * t204 + t203 - (-t136 * t154 + t137 * t195) * t125;
t108 = t184 * t218 + (-qJD(1) * t169 + 0.2e1 * t168 * t189) * t138;
t104 = -0.2e1 * t213 + 0.2e1 * (-t111 * t116 * t128 + (-t116 * t211 - t128 * t213) * t132) * t132;
t1 = [t174 * t194 + (-t146 * t155 * t192 - qJD(4) * t169) * t138, 0, 0, t108, 0, 0; (t120 * t176 + (-t120 * t190 + (-qJD(1) * t110 - t107) * t209) * t114) * t155 + (t121 * t176 * t110 + (((-t113 * t171 + t156 * t184 + t215 * t190) * t136 + (t174 * t199 + t210 + (-t210 + (t216 - t219) * t189) * t138) * t137) * t183 + (-t121 * t190 + t156 * t185) * t110 + (t120 + ((-t149 + t152) * t137 * t181 + t215 * t182) * t121) * t192) * t114) * t157, 0, 0, 0.2e1 * (-t109 * t209 - t120 * t154) * t157 * t214 + ((-t120 * t193 + (-qJD(4) * t109 - t107) * t157 * t121) * t154 + (t120 * t187 + (t108 * t137 * t155 + t217 * t136 - (-qJD(4) * t136 - t113 * t204 + t137 * t191) * t125) * t183 + (-t121 * t193 + t157 * t185) * t109 + ((-t108 - t191) * t136 + ((t125 * t155 - 0.1e1) * qJD(4) + (t125 - t155) * t113) * t137) * t121 * t198) * t156) * t114, 0, 0; (-t127 * t130 + t131 * t206) * t186 + (t131 * t175 - t173 * t127 * t196 + t165 * t207 + (-t173 * t132 * t197 + t131 * t111 - t130 * t112 - t165 * t205) * t128) * t116, 0, 0, t156 * t166 * t186 + (t166 * t190 + (t167 * t193 + ((t127 * t145 + t175) * t144 + (t111 * t144 + (t132 * t145 - t112) * t143) * t128) * t157) * t156) * t116, t104, t104;];
JaD_rot  = t1;
