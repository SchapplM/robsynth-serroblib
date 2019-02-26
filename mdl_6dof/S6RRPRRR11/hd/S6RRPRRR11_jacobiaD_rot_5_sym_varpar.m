% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR11_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:58
% EndTime: 2019-02-26 21:59:59
% DurationCPUTime: 0.68s
% Computational Cost: add. (1381->94), mult. (2734->208), div. (498->12), fcn. (3199->9), ass. (0->96)
t152 = sin(qJ(2));
t145 = 0.1e1 / t152 ^ 2;
t154 = cos(qJ(2));
t149 = t154 ^ 2;
t200 = t145 * t149;
t218 = t154 * t200;
t153 = sin(qJ(1));
t175 = 0.1e1 + t200;
t217 = t153 * t175;
t143 = qJD(4) + qJD(5);
t171 = qJD(1) * t152 + t143;
t155 = cos(qJ(1));
t187 = qJD(2) * t155;
t216 = t171 * t153 - t154 * t187;
t188 = qJD(2) * t154;
t215 = t153 * t188 + t171 * t155;
t194 = t153 * t154;
t136 = atan2(-t194, t152);
t135 = cos(t136);
t134 = sin(t136);
t180 = t134 * t194;
t123 = t135 * t152 - t180;
t120 = 0.1e1 / t123;
t151 = qJ(4) + qJ(5);
t141 = sin(t151);
t142 = cos(t151);
t195 = t153 * t142;
t196 = t152 * t155;
t131 = t141 * t196 + t195;
t127 = 0.1e1 / t131;
t144 = 0.1e1 / t152;
t121 = 0.1e1 / t123 ^ 2;
t128 = 0.1e1 / t131 ^ 2;
t147 = t153 ^ 2;
t140 = t147 * t200 + 0.1e1;
t137 = 0.1e1 / t140;
t214 = t137 - 0.1e1;
t172 = t143 * t152 + qJD(1);
t166 = t172 * t155;
t111 = t141 * t166 + t216 * t142;
t130 = t141 * t153 - t142 * t196;
t126 = t130 ^ 2;
t119 = t126 * t128 + 0.1e1;
t203 = t128 * t130;
t112 = -t216 * t141 + t142 * t166;
t211 = t112 * t127 * t128;
t213 = (t111 * t203 - t126 * t211) / t119 ^ 2;
t191 = qJD(1) * t155;
t178 = t154 * t191;
t189 = qJD(2) * t153;
t113 = ((t152 * t189 - t178) * t144 + t189 * t200) * t137;
t201 = t135 * t154;
t107 = (-t113 * t153 + qJD(2)) * t201 + (-t178 + (-t113 + t189) * t152) * t134;
t212 = t107 * t120 * t121;
t210 = t113 * t134;
t209 = t113 * t154;
t208 = t121 * t154;
t207 = t121 * t155;
t164 = qJD(2) * (-t154 - t218) * t144;
t198 = t149 * t153;
t169 = t191 * t198;
t206 = (t145 * t169 + t147 * t164) / t140 ^ 2;
t125 = t137 * t217;
t205 = t125 * t153;
t204 = t127 * t142;
t202 = t130 * t141;
t150 = t155 ^ 2;
t199 = t149 * t150;
t197 = t152 * t153;
t193 = qJD(1) * t153;
t192 = qJD(1) * t154;
t190 = qJD(2) * t152;
t116 = t121 * t199 + 0.1e1;
t186 = 0.2e1 * (-t199 * t212 + (-t150 * t152 * t188 - t169) * t121) / t116 ^ 2;
t185 = 0.2e1 * t213;
t184 = 0.2e1 * t212;
t183 = -0.2e1 * t206;
t182 = t154 * t207;
t181 = t154 * t206;
t179 = t144 * t198;
t174 = t154 * t186;
t173 = 0.2e1 * t130 * t211;
t170 = t137 * t179;
t168 = t175 * t155;
t167 = t172 * t153;
t165 = t128 * t202 + t204;
t163 = t165 * t155;
t133 = -t141 * t197 + t142 * t155;
t132 = t141 * t155 + t152 * t195;
t117 = 0.1e1 / t119;
t114 = 0.1e1 / t116;
t110 = (t214 * t154 * t134 + t135 * t170) * t155;
t109 = t134 * t197 + t201 + (-t134 * t152 - t135 * t194) * t125;
t108 = t183 * t217 + (qJD(1) * t168 + 0.2e1 * t153 * t164) * t137;
t104 = -0.2e1 * t213 + 0.2e1 * (t111 * t117 * t128 + (-t117 * t211 - t128 * t213) * t130) * t130;
t1 = [0.2e1 * t144 * t155 * t181 + (t144 * t153 * t192 + qJD(2) * t168) * t137, t108, 0, 0, 0, 0; (t120 * t174 + (t120 * t190 + (qJD(1) * t110 + t107) * t208) * t114) * t153 + (t121 * t174 * t110 + (-((-t113 * t170 - t214 * t190 - 0.2e1 * t181) * t134 + (t179 * t183 - t209 + (t209 + (-0.2e1 * t154 - t218) * t189) * t137) * t135) * t182 + (t121 * t190 + t154 * t184) * t110 + (-t120 + ((t147 - t150) * t149 * t144 * t137 * t135 + t214 * t180) * t121) * t192) * t114) * t155 (t109 * t208 + t120 * t152) * t155 * t186 + ((t120 * t193 + (qJD(2) * t109 + t107) * t207) * t152 + (-t120 * t187 - (-t108 * t135 * t153 + t134 * t189 + t205 * t210 - t210 + (-qJD(2) * t134 - t135 * t191) * t125) * t182 + (t121 * t193 + t155 * t184) * t109 - ((-t108 + t191) * t134 + ((-0.1e1 + t205) * qJD(2) + (-t125 + t153) * t113) * t135) * t121 * t196) * t154) * t114, 0, 0, 0, 0; (-t127 * t132 + t133 * t203) * t185 + (t133 * t173 - t127 * t141 * t167 + t215 * t204 + (t130 * t142 * t167 - t133 * t111 - t132 * t112 + t215 * t202) * t128) * t117, t154 * t163 * t185 + (t163 * t190 + (t165 * t193 + ((t127 * t143 + t173) * t141 + (-t111 * t141 + (-t130 * t143 + t112) * t142) * t128) * t155) * t154) * t117, 0, t104, t104, 0;];
JaD_rot  = t1;
