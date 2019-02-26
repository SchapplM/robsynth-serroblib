% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:22
% EndTime: 2019-02-26 21:35:23
% DurationCPUTime: 1.05s
% Computational Cost: add. (3877->124), mult. (6168->272), div. (1114->15), fcn. (7752->9), ass. (0->115)
t162 = qJ(2) + pkin(9);
t161 = cos(t162);
t167 = sin(qJ(4));
t241 = sin(qJ(1));
t201 = t241 * t167;
t168 = cos(qJ(4));
t169 = cos(qJ(1));
t221 = t169 * t168;
t145 = t161 * t221 + t201;
t139 = 0.1e1 / t145 ^ 2;
t160 = sin(t162);
t156 = t160 ^ 2;
t166 = t169 ^ 2;
t229 = t156 * t166;
t206 = t139 * t229;
t134 = 0.1e1 + t206;
t193 = qJD(1) * t241;
t218 = qJD(2) * t169;
t197 = t160 * t218;
t179 = t161 * t193 + t197;
t192 = t241 * qJD(4);
t222 = t169 * t167;
t124 = (-qJD(4) * t161 + qJD(1)) * t222 + (t192 - t179) * t168;
t138 = 0.1e1 / t145;
t236 = t124 * t138 * t139;
t187 = t229 * t236;
t198 = qJD(2) * t160 * t166;
t244 = (-t187 + (-t156 * t169 * t193 + t161 * t198) * t139) / t134 ^ 2;
t224 = t160 * t169;
t141 = t161 * t201 + t221;
t184 = t167 * t192;
t215 = qJD(4) * t169;
t195 = t168 * t215;
t123 = t141 * qJD(1) - t161 * t195 + t167 * t197 - t184;
t200 = t241 * t168;
t144 = t161 * t222 - t200;
t157 = 0.1e1 / t160;
t163 = 0.1e1 / t167;
t164 = 0.1e1 / t167 ^ 2;
t216 = qJD(4) * t168;
t196 = t164 * t216;
t158 = 0.1e1 / t160 ^ 2;
t219 = qJD(2) * t161;
t199 = t158 * t219;
t228 = t157 * t163;
t243 = (t157 * t196 + t163 * t199) * t144 + t123 * t228;
t225 = t160 * t167;
t133 = atan2(-t141, t225);
t128 = cos(t133);
t127 = sin(t133);
t235 = t127 * t141;
t122 = t128 * t225 - t235;
t119 = 0.1e1 / t122;
t120 = 0.1e1 / t122 ^ 2;
t242 = 0.2e1 * t144;
t136 = t141 ^ 2;
t226 = t158 * t164;
t135 = t136 * t226 + 0.1e1;
t131 = 0.1e1 / t135;
t180 = t160 * t216 + t167 * t219;
t204 = t141 * t226;
t202 = t241 * t160;
t185 = qJD(2) * t202;
t186 = t168 * t193;
t220 = qJD(1) * t169;
t125 = t168 * t192 * t161 - t186 + (t220 * t161 - t185 - t215) * t167;
t207 = t125 * t228;
t111 = (t180 * t204 - t207) * t131;
t177 = -t111 * t141 + t180;
t107 = (-t111 * t225 - t125) * t127 + t177 * t128;
t121 = t119 * t120;
t240 = t107 * t121;
t159 = t157 / t156;
t165 = t163 * t164;
t239 = (t125 * t204 + (-t158 * t165 * t216 - t159 * t164 * t219) * t136) / t135 ^ 2;
t238 = t120 * t144;
t237 = t123 * t120;
t234 = t127 * t144;
t233 = t127 * t160;
t232 = t128 * t141;
t231 = t128 * t144;
t230 = t128 * t161;
t227 = t158 * t161;
t223 = t164 * t168;
t217 = qJD(4) * t167;
t137 = t144 ^ 2;
t117 = t120 * t137 + 0.1e1;
t214 = 0.2e1 * (-t137 * t240 - t144 * t237) / t117 ^ 2;
t213 = 0.2e1 * t244;
t212 = -0.2e1 * t239;
t211 = t121 * t242;
t210 = t157 * t239;
t209 = t120 * t234;
t205 = t141 * t228;
t203 = t163 * t227;
t182 = t141 * t203 + t241;
t118 = t182 * t131;
t194 = t241 - t118;
t191 = t119 * t214;
t190 = t120 * t214;
t189 = t224 * t242;
t188 = t163 * t210;
t143 = t161 * t200 - t222;
t183 = t141 * t223 - t143 * t163;
t181 = t139 * t143 * t169 - t241 * t138;
t129 = 0.1e1 / t134;
t126 = t145 * qJD(1) - t161 * t184 - t168 * t185 - t195;
t115 = 0.1e1 / t117;
t114 = t183 * t157 * t131;
t110 = (-t127 + (t128 * t205 + t127) * t131) * t144;
t109 = -t118 * t232 + (t194 * t233 + t230) * t167;
t108 = t128 * t160 * t168 - t127 * t143 + (-t127 * t225 - t232) * t114;
t106 = t182 * t212 + (t125 * t203 + t220 + (-t196 * t227 + (-0.2e1 * t159 * t161 ^ 2 - t157) * t163 * qJD(2)) * t141) * t131;
t104 = -0.2e1 * t183 * t210 + (-t183 * t199 + (t125 * t223 - t126 * t163 + (t143 * t223 + (-0.2e1 * t165 * t168 ^ 2 - t163) * t141) * qJD(4)) * t157) * t131;
t1 = [t243 * t131 + t188 * t242, t106, 0, t104, 0, 0; t141 * t191 + (-t125 * t119 + (t107 * t141 + t110 * t123) * t120) * t115 + (t110 * t190 + (0.2e1 * t110 * t240 + (t123 * t131 - t123 - (-t111 * t131 * t205 + t212) * t144) * t120 * t127 + (-(-0.2e1 * t141 * t188 - t111) * t238 + (-(t111 + t207) * t144 + t243 * t141) * t120 * t131) * t128) * t115) * t144, t109 * t144 * t190 + (-(-t106 * t232 + (t111 * t235 - t125 * t128) * t118) * t238 + (t107 * t211 + t237) * t109 + (-t119 * t224 - (-t118 * t233 + t127 * t202 + t230) * t238) * t216) * t115 + (t191 * t224 + ((-t119 * t218 - (t194 * qJD(2) - t111) * t209) * t161 + (t119 * t193 + (t169 * t107 - (-t106 + t220) * t234 - (t194 * t111 - qJD(2)) * t231) * t120) * t160) * t115) * t167, 0 (t108 * t238 - t119 * t145) * t214 + (t108 * t237 + t124 * t119 + (t108 * t211 - t120 * t145) * t107 - (t168 * t219 - t160 * t217 - t104 * t141 - t114 * t125 + (-t114 * t225 - t143) * t111) * t120 * t231 - (-t126 + (-t104 * t167 - t111 * t168) * t160 - t177 * t114) * t209) * t115, 0, 0; t181 * t160 * t213 + (-t181 * t219 + ((qJD(1) * t138 + 0.2e1 * t143 * t236) * t169 + (-t241 * t124 - t126 * t169 + t143 * t193) * t139) * t160) * t129 (t138 * t161 * t169 + t168 * t206) * t213 + (0.2e1 * t168 * t187 + t179 * t138 + ((t124 * t169 - 0.2e1 * t168 * t198) * t161 + (t166 * t217 + 0.2e1 * t169 * t186) * t156) * t139) * t129, 0, t139 * t189 * t244 + (t189 * t236 + (t123 * t224 + (t160 * t193 - t161 * t218) * t144) * t139) * t129, 0, 0;];
JaD_rot  = t1;
