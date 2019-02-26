% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:48
% EndTime: 2019-02-26 22:18:49
% DurationCPUTime: 0.76s
% Computational Cost: add. (2009->96), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
t154 = sin(qJ(2));
t148 = t154 ^ 2;
t156 = cos(qJ(2));
t151 = 0.1e1 / t156 ^ 2;
t201 = t148 * t151;
t155 = sin(qJ(1));
t219 = 0.2e1 * t155;
t218 = t154 * t201;
t145 = qJ(3) + pkin(11) + qJ(5);
t144 = cos(t145);
t157 = cos(qJ(1));
t193 = t156 * t157;
t143 = sin(t145);
t197 = t155 * t143;
t133 = t144 * t193 + t197;
t195 = t155 * t154;
t138 = atan2(-t195, -t156);
t136 = cos(t138);
t135 = sin(t138);
t182 = t135 * t195;
t123 = -t136 * t156 - t182;
t120 = 0.1e1 / t123;
t127 = 0.1e1 / t133;
t150 = 0.1e1 / t156;
t121 = 0.1e1 / t123 ^ 2;
t128 = 0.1e1 / t133 ^ 2;
t217 = -0.2e1 * t154;
t149 = t155 ^ 2;
t141 = t149 * t201 + 0.1e1;
t139 = 0.1e1 / t141;
t216 = t139 - 0.1e1;
t146 = qJD(3) + qJD(5);
t194 = t155 * t156;
t167 = t143 * t194 + t144 * t157;
t188 = qJD(2) * t157;
t178 = t154 * t188;
t111 = t167 * qJD(1) - t133 * t146 + t143 * t178;
t196 = t155 * t144;
t132 = t143 * t193 - t196;
t126 = t132 ^ 2;
t116 = t126 * t128 + 0.1e1;
t206 = t128 * t132;
t172 = -qJD(1) * t156 + t146;
t173 = t146 * t156 - qJD(1);
t203 = t143 * t157;
t112 = -t173 * t203 + (t172 * t155 - t178) * t144;
t213 = t112 * t127 * t128;
t215 = (-t111 * t206 - t126 * t213) / t116 ^ 2;
t191 = qJD(1) * t157;
t179 = t154 * t191;
t189 = qJD(2) * t156;
t190 = qJD(2) * t155;
t113 = (-(-t155 * t189 - t179) * t150 + t190 * t201) * t139;
t204 = t136 * t154;
t107 = (-t113 * t155 + qJD(2)) * t204 + (-t179 + (t113 - t190) * t156) * t135;
t214 = t107 * t120 * t121;
t212 = t113 * t135;
t211 = t113 * t154;
t210 = t121 * t154;
t199 = t150 * t154;
t166 = qJD(2) * (t150 * t218 + t199);
t170 = t148 * t155 * t191;
t209 = (t149 * t166 + t151 * t170) / t141 ^ 2;
t177 = 0.1e1 + t201;
t125 = t177 * t155 * t139;
t208 = t125 * t155;
t207 = t127 * t143;
t205 = t132 * t144;
t202 = t148 * t150;
t153 = t157 ^ 2;
t200 = t148 * t153;
t198 = t154 * t157;
t192 = qJD(1) * t155;
t119 = t121 * t200 + 0.1e1;
t187 = 0.2e1 * (-t200 * t214 + (t153 * t154 * t189 - t170) * t121) / t119 ^ 2;
t186 = -0.2e1 * t215;
t185 = 0.2e1 * t214;
t184 = t132 * t213;
t183 = t121 * t198;
t181 = t139 * t202;
t176 = t154 * t187;
t175 = t209 * t217;
t174 = t209 * t219;
t171 = t155 * t181;
t169 = t177 * t157;
t168 = t128 * t205 - t207;
t165 = t154 * t190 + t172 * t157;
t131 = -t144 * t194 + t203;
t117 = 0.1e1 / t119;
t114 = 0.1e1 / t116;
t110 = (t216 * t154 * t135 - t136 * t171) * t157;
t109 = -t135 * t194 + t204 + (t135 * t156 - t136 * t195) * t125;
t108 = -t177 * t174 + (qJD(1) * t169 + t166 * t219) * t139;
t104 = t186 + 0.2e1 * (-t111 * t114 * t128 + (-t114 * t213 - t128 * t215) * t132) * t132;
t1 = [t150 * t157 * t175 + (qJD(2) * t169 - t192 * t199) * t139, t108, 0, 0, 0, 0; (t120 * t176 + (-t120 * t189 + (qJD(1) * t110 + t107) * t210) * t117) * t155 + (t121 * t176 * t110 + (-((t113 * t171 + t216 * t189 + t175) * t135 + (t174 * t202 - t211 + (t211 + (t217 - t218) * t190) * t139) * t136) * t183 + (-t121 * t189 + t154 * t185) * t110 + (-t120 + ((-t149 + t153) * t136 * t181 + t216 * t182) * t121) * t154 * qJD(1)) * t117) * t157 (t109 * t210 - t120 * t156) * t157 * t187 + ((-t120 * t192 + (-qJD(2) * t109 - t107) * t157 * t121) * t156 + (-t120 * t188 - (-t108 * t136 * t155 + t135 * t190 + t208 * t212 - t212 + (-qJD(2) * t135 - t136 * t191) * t125) * t183 + (t121 * t192 + t157 * t185) * t109 - ((t108 - t191) * t135 + ((0.1e1 - t208) * qJD(2) + (t125 - t155) * t113) * t136) * t121 * t193) * t154) * t117, 0, 0, 0, 0; 0.2e1 * (t127 * t167 + t131 * t206) * t215 + (0.2e1 * t131 * t184 - t173 * t127 * t196 + t165 * t207 + (-t173 * t132 * t197 + t131 * t111 + t112 * t167 - t165 * t205) * t128) * t114, t168 * t186 * t198 + (t168 * t156 * t188 + (-t168 * t192 + ((-t127 * t146 - 0.2e1 * t184) * t144 + (-t111 * t144 + (-t132 * t146 + t112) * t143) * t128) * t157) * t154) * t114, t104, 0, t104, 0;];
JaD_rot  = t1;
