% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:57:35
% EndTime: 2019-02-26 20:57:36
% DurationCPUTime: 1.00s
% Computational Cost: add. (3939->118), mult. (5988->261), div. (1156->14), fcn. (7606->9), ass. (0->107)
t146 = qJ(1) + pkin(9);
t144 = sin(t146);
t145 = cos(t146);
t156 = cos(qJ(4));
t154 = sin(qJ(4));
t157 = cos(qJ(3));
t208 = t154 * t157;
t126 = t144 * t208 + t145 * t156;
t155 = sin(qJ(3));
t207 = t155 * t154;
t121 = atan2(-t126, t207);
t117 = sin(t121);
t118 = cos(t121);
t123 = t126 ^ 2;
t148 = 0.1e1 / t154 ^ 2;
t151 = 0.1e1 / t155 ^ 2;
t211 = t148 * t151;
t122 = t123 * t211 + 0.1e1;
t119 = 0.1e1 / t122;
t147 = 0.1e1 / t154;
t209 = t151 * t157;
t187 = t147 * t209;
t171 = t126 * t187 + t144;
t105 = t171 * t119;
t205 = -t105 + t144;
t229 = t205 * t155 * t117 + t118 * t157;
t150 = 0.1e1 / t155;
t152 = t150 * t151;
t227 = qJD(3) * (0.2e1 * t152 * t157 ^ 2 + t150);
t199 = qJD(4) * t156;
t179 = t145 * t199;
t200 = qJD(4) * t154;
t180 = t144 * t200;
t202 = qJD(3) * t155;
t183 = t154 * t202;
t110 = qJD(1) * t126 + t145 * t183 - t157 * t179 - t180;
t214 = t144 * t156;
t129 = t145 * t208 - t214;
t201 = qJD(3) * t157;
t185 = t151 * t201;
t212 = t147 * t150;
t226 = t129 * (t148 * t150 * t199 + t147 * t185) + t110 * t212;
t222 = t117 * t126;
t109 = t118 * t207 - t222;
t106 = 0.1e1 / t109;
t141 = 0.1e1 / t145;
t107 = 0.1e1 / t109 ^ 2;
t142 = 0.1e1 / t145 ^ 2;
t168 = t154 * t201 + t155 * t199;
t189 = t126 * t211;
t203 = qJD(1) * t145;
t204 = qJD(1) * t144;
t112 = -t144 * t183 - t145 * t200 - t156 * t204 + (t144 * t199 + t154 * t203) * t157;
t191 = t112 * t212;
t98 = (t168 * t189 - t191) * t119;
t166 = -t126 * t98 + t168;
t172 = -t98 * t207 - t112;
t94 = t117 * t172 + t118 * t166;
t225 = t106 * t107 * t94;
t149 = t147 * t148;
t181 = t151 * t199;
t184 = t152 * t201;
t224 = (t112 * t189 + (-t148 * t184 - t149 * t181) * t123) / t122 ^ 2;
t223 = t107 * t129;
t221 = t117 * t129;
t219 = t118 * t126;
t218 = t118 * t129;
t216 = t142 * t144;
t215 = t142 * t151;
t213 = t145 * t106;
t210 = t148 * t156;
t206 = t156 * t157;
t124 = t129 ^ 2;
t104 = t107 * t124 + 0.1e1;
t198 = 0.2e1 / t104 ^ 2 * (-t110 * t223 - t124 * t225);
t197 = 0.2e1 * t225;
t182 = t156 * t202;
t111 = (-qJD(1) * t157 + qJD(4)) * t214 + (-t182 + (-qJD(4) * t157 + qJD(1)) * t154) * t145;
t130 = t144 * t154 + t145 * t206;
t125 = t130 ^ 2;
t116 = t125 * t215 + 0.1e1;
t143 = t141 * t142;
t186 = t151 * t204;
t196 = 0.2e1 / t116 ^ 2 * (t130 * t111 * t215 + (-t142 * t184 + t143 * t186) * t125);
t195 = -0.2e1 * t224;
t194 = t150 * t224;
t193 = t107 * t221;
t190 = t126 * t212;
t188 = t141 * t209;
t178 = t106 * t198;
t177 = t107 * t198;
t176 = t129 * t197;
t175 = t150 * t196;
t173 = t147 * t194;
t128 = t144 * t206 - t145 * t154;
t170 = t126 * t210 - t128 * t147;
t169 = t128 * t141 - t130 * t216;
t114 = 0.1e1 / t116;
t113 = qJD(1) * t130 - t144 * t182 - t157 * t180 - t179;
t102 = 0.1e1 / t104;
t101 = t170 * t150 * t119;
t97 = (-t117 + (t118 * t190 + t117) * t119) * t129;
t96 = -t105 * t219 + t229 * t154;
t95 = t118 * t155 * t156 - t117 * t128 + (-t117 * t207 - t219) * t101;
t93 = t171 * t195 + (t112 * t187 + t203 + (-t148 * t157 * t181 - t147 * t227) * t126) * t119;
t91 = -0.2e1 * t170 * t194 + (-t170 * t185 + (t112 * t210 - t113 * t147 + (t128 * t210 + (-0.2e1 * t149 * t156 ^ 2 - t147) * t126) * qJD(4)) * t150) * t119;
t1 = [t226 * t119 + 0.2e1 * t129 * t173, 0, t93, t91, 0, 0; t126 * t178 + (-t112 * t106 + (t110 * t97 + t126 * t94) * t107) * t102 + (t97 * t177 + (t97 * t197 + (t110 * t119 - t110 - (-t119 * t98 * t190 + t195) * t129) * t107 * t117 + (-(-0.2e1 * t126 * t173 - t98) * t223 + (-(t98 + t191) * t129 + t226 * t126) * t107 * t119) * t118) * t102) * t129, 0, t96 * t129 * t177 + (t96 * t176 + (-(-t93 * t219 + (-t112 * t118 + t222 * t98) * t105) * t129 + t96 * t110) * t107 + (-t155 * t213 - t229 * t223) * t199) * t102 + (t145 * t178 * t155 + ((-qJD(3) * t213 - (t205 * qJD(3) - t98) * t193) * t157 + (t106 * t204 + (t145 * t94 - (-t93 + t203) * t221 - (t205 * t98 - qJD(3)) * t218) * t107) * t155) * t102) * t154 (-t106 * t130 + t95 * t223) * t198 + (t95 * t176 + t111 * t106 - (-t113 + (-t154 * t91 - t156 * t98) * t155 - t166 * t101) * t193 + (t95 * t110 - t130 * t94 - (t101 * t172 - t126 * t91 - t128 * t98 - t155 * t200 + t156 * t201) * t218) * t107) * t102, 0, 0; t169 * t175 + (t169 * t185 + (t111 * t216 - t113 * t141 + (-t128 * t216 + (0.2e1 * t143 * t144 ^ 2 + t141) * t130) * qJD(1)) * t150) * t114, 0 (t130 * t188 + t156) * t196 + (-t111 * t188 + t200 + (-t142 * t157 * t186 + t141 * t227) * t130) * t114, t129 * t141 * t175 + (t110 * t141 * t150 + (-t142 * t150 * t204 + t141 * t185) * t129) * t114, 0, 0;];
JaD_rot  = t1;
