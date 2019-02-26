% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP3
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
% Datum: 2019-02-26 21:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:56
% EndTime: 2019-02-26 21:35:57
% DurationCPUTime: 0.96s
% Computational Cost: add. (4725->120), mult. (5988->260), div. (1156->14), fcn. (7606->9), ass. (0->113)
t154 = sin(qJ(2));
t157 = cos(qJ(1));
t232 = t154 * t157;
t147 = 0.1e1 / t154;
t148 = 0.1e1 / t154 ^ 2;
t149 = t147 * t148;
t156 = cos(qJ(2));
t231 = qJD(2) * (0.2e1 * t149 * t156 ^ 2 + t147);
t146 = pkin(9) + qJ(4);
t144 = sin(t146);
t145 = cos(t146);
t155 = sin(qJ(1));
t198 = qJD(4) * t156;
t180 = t144 * t198;
t202 = qJD(2) * t157;
t181 = t154 * t202;
t205 = qJD(1) * t156;
t185 = t155 * t205;
t200 = qJD(4) * t145;
t204 = qJD(1) * t157;
t115 = t157 * t180 - t155 * t200 - t144 * t204 + (t181 + t185) * t145;
t216 = t145 * t157;
t134 = t155 * t144 + t156 * t216;
t141 = 0.1e1 / t145;
t142 = 0.1e1 / t145 ^ 2;
t203 = qJD(2) * t156;
t184 = t148 * t203;
t201 = qJD(4) * t144;
t219 = t141 * t147;
t230 = (-t142 * t147 * t201 + t141 * t184) * t134 + t115 * t219;
t208 = t157 * t144;
t209 = t155 * t156;
t131 = t145 * t209 - t208;
t212 = t154 * t145;
t122 = atan2(-t131, t212);
t119 = cos(t122);
t118 = sin(t122);
t225 = t118 * t131;
t113 = t119 * t212 - t225;
t110 = 0.1e1 / t113;
t151 = 0.1e1 / t157;
t111 = 0.1e1 / t113 ^ 2;
t152 = 0.1e1 / t157 ^ 2;
t127 = t131 ^ 2;
t217 = t142 * t148;
t123 = t127 * t217 + 0.1e1;
t120 = 0.1e1 / t123;
t199 = qJD(4) * t154;
t168 = -t144 * t199 + t145 * t203;
t188 = t131 * t217;
t211 = t154 * t155;
t182 = qJD(2) * t211;
t117 = t134 * qJD(1) - t145 * t182 - t155 * t180 - t157 * t200;
t190 = t117 * t219;
t102 = (t168 * t188 - t190) * t120;
t166 = -t102 * t131 + t168;
t98 = (-t102 * t212 - t117) * t118 + t166 * t119;
t229 = t110 * t111 * t98;
t143 = t141 * t142;
t183 = t149 * t203;
t228 = (t117 * t188 + (t143 * t148 * t201 - t142 * t183) * t127) / t123 ^ 2;
t227 = t111 * t134;
t226 = t115 * t111;
t224 = t118 * t134;
t223 = t118 * t154;
t222 = t119 * t131;
t221 = t119 * t134;
t220 = t119 * t156;
t218 = t142 * t144;
t215 = t148 * t152;
t214 = t148 * t156;
t213 = t152 * t155;
t210 = t155 * t145;
t187 = t141 * t214;
t171 = t131 * t187 + t155;
t109 = t171 * t120;
t207 = -t109 + t155;
t206 = qJD(1) * t155;
t129 = t134 ^ 2;
t108 = t111 * t129 + 0.1e1;
t197 = 0.2e1 / t108 ^ 2 * (-t129 * t229 - t134 * t226);
t196 = 0.2e1 * t229;
t195 = -0.2e1 * t228;
t173 = -qJD(4) + t205;
t174 = -qJD(1) + t198;
t114 = -t174 * t216 + (t173 * t155 + t181) * t144;
t133 = -t156 * t208 + t210;
t128 = t133 ^ 2;
t126 = t128 * t215 + 0.1e1;
t153 = t151 * t152;
t194 = 0.2e1 * (t133 * t114 * t215 + (t148 * t153 * t206 - t152 * t183) * t128) / t126 ^ 2;
t193 = t147 * t228;
t192 = t111 * t224;
t189 = t131 * t219;
t186 = t151 * t214;
t179 = t110 * t197;
t178 = t111 * t197;
t177 = t134 * t196;
t176 = t147 * t194;
t172 = t141 * t193;
t130 = t144 * t209 + t216;
t170 = -t130 * t141 + t131 * t218;
t169 = -t130 * t151 - t133 * t213;
t124 = 0.1e1 / t126;
t116 = t174 * t210 + (t173 * t157 - t182) * t144;
t106 = 0.1e1 / t108;
t105 = t170 * t147 * t120;
t101 = (-t118 + (t119 * t189 + t118) * t120) * t134;
t100 = -t109 * t222 + (t207 * t223 + t220) * t145;
t99 = -t119 * t154 * t144 + t118 * t130 - (-t118 * t212 - t222) * t105;
t97 = t171 * t195 + (t117 * t187 + t204 + (-t141 * t231 + t180 * t217) * t131) * t120;
t95 = 0.2e1 * t170 * t193 + (t170 * t184 + (-t117 * t218 + t116 * t141 + (t130 * t218 + (-0.2e1 * t143 * t144 ^ 2 - t141) * t131) * qJD(4)) * t147) * t120;
t1 = [t120 * t230 + 0.2e1 * t134 * t172, t97, 0, t95, 0, 0; t131 * t179 + (-t117 * t110 + (t101 * t115 + t131 * t98) * t111) * t106 + (t101 * t178 + (t101 * t196 + (t115 * t120 - t115 - (-t102 * t120 * t189 + t195) * t134) * t111 * t118 + (-(-0.2e1 * t131 * t172 - t102) * t227 + (-(t102 + t190) * t134 + t230 * t131) * t111 * t120) * t119) * t106) * t134, t100 * t134 * t178 + (-(-t97 * t222 + (t102 * t225 - t117 * t119) * t109) * t227 + (t177 + t226) * t100 + (t110 * t232 - (t109 * t223 - t118 * t211 - t220) * t227) * t201) * t106 + (t179 * t232 + ((-t110 * t202 - (t207 * qJD(2) - t102) * t192) * t156 + (t110 * t206 + (t157 * t98 - (-t97 + t204) * t224 - (t207 * t102 - qJD(2)) * t221) * t111) * t154) * t106) * t145, 0 (-t110 * t133 + t99 * t227) * t197 + (t99 * t177 + t114 * t110 - (t116 + (t102 * t144 - t145 * t95) * t154 + t166 * t105) * t192 + (t99 * t115 - t133 * t98 - (-t144 * t203 - t145 * t199 + t105 * t117 - t131 * t95 + (t105 * t212 + t130) * t102) * t221) * t111) * t106, 0, 0; t169 * t176 + (t169 * t184 + (t114 * t213 + t116 * t151 + (t130 * t213 + (0.2e1 * t153 * t155 ^ 2 + t151) * t133) * qJD(1)) * t147) * t124 (t133 * t186 - t144) * t194 + (-t114 * t186 + t200 + (t151 * t231 - t185 * t215) * t133) * t124, 0, t134 * t151 * t176 + (t115 * t147 * t151 + (-t147 * t152 * t206 + t151 * t184) * t134) * t124, 0, 0;];
JaD_rot  = t1;
