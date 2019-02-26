% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:47:36
% EndTime: 2019-02-26 20:47:37
% DurationCPUTime: 0.95s
% Computational Cost: add. (3877->122), mult. (6168->268), div. (1114->15), fcn. (7752->9), ass. (0->115)
t156 = qJ(3) + pkin(9);
t155 = cos(t156);
t162 = sin(qJ(1));
t217 = t155 * t162;
t154 = sin(t156);
t161 = sin(qJ(5));
t236 = cos(qJ(1));
t196 = t236 * t161;
t163 = cos(qJ(5));
t214 = t162 * t163;
t142 = t154 * t196 + t214;
t218 = t155 * t161;
t132 = atan2(t142, t218);
t127 = cos(t132);
t126 = sin(t132);
t227 = t126 * t142;
t121 = t127 * t218 + t227;
t118 = 0.1e1 / t121;
t141 = t154 * t214 + t196;
t136 = 0.1e1 / t141;
t151 = 0.1e1 / t155;
t157 = 0.1e1 / t161;
t119 = 0.1e1 / t121 ^ 2;
t137 = 0.1e1 / t141 ^ 2;
t158 = 0.1e1 / t161 ^ 2;
t195 = t236 * t163;
t215 = t162 * t161;
t140 = t154 * t215 - t195;
t135 = t140 ^ 2;
t116 = t119 * t135 + 0.1e1;
t189 = qJD(1) * t236;
t179 = t154 * t189;
t173 = t236 * qJD(5) + t179;
t184 = qJD(5) * t154 + qJD(1);
t211 = qJD(3) * t162;
t193 = t155 * t211;
t124 = t184 * t214 + (t173 + t193) * t161;
t230 = t124 * t119;
t139 = t142 ^ 2;
t152 = 0.1e1 / t155 ^ 2;
t219 = t152 * t158;
t134 = t139 * t219 + 0.1e1;
t130 = 0.1e1 / t134;
t208 = qJD(5) * t163;
t212 = qJD(3) * t154;
t174 = t155 * t208 - t161 * t212;
t199 = t142 * t219;
t197 = t236 * t155;
t180 = qJD(3) * t197;
t181 = t163 * t189;
t182 = t154 * t195;
t122 = -qJD(5) * t182 - t161 * t180 - t181 + (qJD(1) * t154 + qJD(5)) * t215;
t221 = t151 * t157;
t202 = t122 * t221;
t110 = (-t174 * t199 - t202) * t130;
t172 = -t110 * t142 - t174;
t106 = (-t110 * t218 - t122) * t126 - t172 * t127;
t120 = t118 * t119;
t234 = t106 * t120;
t235 = (-t135 * t234 + t140 * t230) / t116 ^ 2;
t150 = t155 ^ 2;
t160 = t162 ^ 2;
t222 = t150 * t160;
t201 = t137 * t222;
t133 = 0.1e1 + t201;
t210 = qJD(3) * t163;
t192 = t155 * t210;
t125 = t173 * t163 + (-t161 * t184 + t192) * t162;
t229 = t125 * t136 * t137;
t183 = t222 * t229;
t233 = (-t183 + (t150 * t162 * t189 - t155 * t160 * t212) * t137) / t133 ^ 2;
t153 = t151 / t150;
t159 = t157 * t158;
t232 = (-t122 * t199 + (-t152 * t159 * t208 + t153 * t158 * t212) * t139) / t134 ^ 2;
t231 = t119 * t140;
t228 = t126 * t140;
t226 = t126 * t155;
t225 = t127 * t140;
t224 = t127 * t142;
t223 = t127 * t154;
t220 = t152 * t154;
t216 = t158 * t163;
t213 = qJD(1) * t162;
t209 = qJD(5) * t161;
t207 = 0.2e1 * t235;
t206 = 0.2e1 * t233;
t205 = -0.2e1 * t232;
t204 = 0.2e1 * t120 * t140;
t203 = t119 * t228;
t200 = t142 * t221;
t198 = t157 * t220;
t194 = t152 * t212;
t191 = t158 * t208;
t176 = t142 * t198 + t236;
t117 = t176 * t130;
t190 = t236 - t117;
t188 = -0.2e1 * t118 * t235;
t187 = t119 * t207;
t186 = 0.2e1 * t151 * t232;
t185 = -0.2e1 * t140 * t217;
t178 = t157 * t186;
t143 = t182 - t215;
t177 = t142 * t216 - t143 * t157;
t175 = t137 * t143 * t162 - t236 * t136;
t171 = t124 * t221 - (t151 * t191 - t157 * t194) * t140;
t128 = 0.1e1 / t133;
t123 = qJD(1) * t141 + qJD(5) * t142 - t163 * t180;
t114 = 0.1e1 / t116;
t113 = t177 * t151 * t130;
t109 = (-t126 + (-t127 * t200 + t126) * t130) * t140;
t108 = t117 * t224 + (t190 * t226 - t223) * t161;
t107 = t127 * t155 * t163 + t126 * t143 - (-t126 * t218 + t224) * t113;
t105 = t176 * t205 + (-t122 * t198 - t213 + (-t191 * t220 + (0.2e1 * t153 * t154 ^ 2 + t151) * t157 * qJD(3)) * t142) * t130;
t103 = t177 * t186 + (-t177 * t194 + (t122 * t216 - t123 * t157 + (-t143 * t216 + (0.2e1 * t159 * t163 ^ 2 + t157) * t142) * qJD(5)) * t151) * t130;
t1 = [-t130 * t171 + t140 * t178, 0, t105, 0, t103, 0; t142 * t188 + (-t122 * t118 + (-t106 * t142 - t109 * t124) * t119) * t114 + (t109 * t187 + (0.2e1 * t109 * t234 + (-t124 * t130 + t124 - (t110 * t130 * t200 + t205) * t140) * t119 * t126 + (-(t142 * t178 - t110) * t231 + (-(t110 + t202) * t140 + t171 * t142) * t119 * t130) * t127) * t114) * t140, 0, t108 * t140 * t187 + (-(t105 * t224 + (-t110 * t227 - t122 * t127) * t117) * t231 + (t106 * t204 - t230) * t108 + (t118 * t217 - (-t117 * t226 + t126 * t197 - t223) * t231) * t208) * t114 + (t188 * t217 + ((-t118 * t211 - (-qJD(3) * t190 + t110) * t203) * t154 + (t118 * t189 + (-t162 * t106 - (-t105 - t213) * t228 - (t110 * t190 - qJD(3)) * t225) * t119) * t155) * t114) * t161, 0 (t107 * t231 - t118 * t141) * t207 + (-t107 * t230 + t125 * t118 + (t107 * t204 - t119 * t141) * t106 - (-t154 * t210 - t155 * t209 + t103 * t142 + t113 * t122 + (t113 * t218 + t143) * t110) * t119 * t225 - (-t123 + (-t103 * t161 - t110 * t163) * t155 - t172 * t113) * t203) * t114, 0; t175 * t155 * t206 + (t175 * t212 + ((-qJD(1) * t136 + 0.2e1 * t143 * t229) * t162 + (t123 * t162 - t236 * t125 - t143 * t189) * t137) * t155) * t128, 0 (t136 * t154 * t162 + t163 * t201) * t206 + (0.2e1 * t163 * t183 + (-t179 - t193) * t136 + ((t125 * t162 + 0.2e1 * t160 * t192) * t154 + (t160 * t209 - 0.2e1 * t162 * t181) * t150) * t137) * t128, 0, t137 * t185 * t233 + (t185 * t229 + (t124 * t217 + (-t154 * t211 + t155 * t189) * t140) * t137) * t128, 0;];
JaD_rot  = t1;
