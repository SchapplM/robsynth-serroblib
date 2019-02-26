% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_rot = S6RRPRPP3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:56
% EndTime: 2019-02-26 21:35:57
% DurationCPUTime: 1.02s
% Computational Cost: add. (4725->120), mult. (5988->261), div. (1156->14), fcn. (7606->9), ass. (0->111)
t159 = sin(qJ(2));
t162 = cos(qJ(1));
t235 = t159 * t162;
t152 = 0.1e1 / t159;
t153 = 0.1e1 / t159 ^ 2;
t154 = t152 * t153;
t161 = cos(qJ(2));
t234 = qJD(2) * (0.2e1 * t154 * t161 ^ 2 + t152);
t151 = pkin(9) + qJ(4);
t149 = sin(t151);
t160 = sin(qJ(1));
t213 = t160 * t161;
t150 = cos(t151);
t219 = t150 * t162;
t131 = t149 * t213 + t219;
t203 = qJD(4) * t162;
t183 = t150 * t203;
t204 = qJD(4) * t160;
t184 = t149 * t204;
t207 = qJD(2) * t162;
t186 = t159 * t207;
t115 = t131 * qJD(1) + t149 * t186 - t161 * t183 - t184;
t212 = t162 * t149;
t134 = -t160 * t150 + t161 * t212;
t146 = 0.1e1 / t149;
t147 = 0.1e1 / t149 ^ 2;
t208 = qJD(2) * t161;
t189 = t153 * t208;
t206 = qJD(4) * t150;
t222 = t146 * t152;
t233 = (t147 * t152 * t206 + t146 * t189) * t134 + t115 * t222;
t215 = t159 * t149;
t123 = atan2(-t131, t215);
t120 = cos(t123);
t119 = sin(t123);
t228 = t119 * t131;
t114 = t120 * t215 - t228;
t111 = 0.1e1 / t114;
t156 = 0.1e1 / t162;
t112 = 0.1e1 / t114 ^ 2;
t157 = 0.1e1 / t162 ^ 2;
t128 = t131 ^ 2;
t220 = t147 * t153;
t124 = t128 * t220 + 0.1e1;
t121 = 0.1e1 / t124;
t205 = qJD(4) * t159;
t173 = t149 * t208 + t150 * t205;
t193 = t131 * t220;
t214 = t159 * t160;
t187 = qJD(2) * t214;
t209 = qJD(1) * t162;
t210 = qJD(1) * t160;
t117 = (t204 * t161 - t210) * t150 + (t209 * t161 - t187 - t203) * t149;
t195 = t117 * t222;
t103 = (t173 * t193 - t195) * t121;
t171 = -t103 * t131 + t173;
t99 = (-t103 * t215 - t117) * t119 + t171 * t120;
t232 = t111 * t112 * t99;
t148 = t146 * t147;
t185 = t153 * t206;
t188 = t154 * t208;
t231 = (t117 * t193 + (-t147 * t188 - t148 * t185) * t128) / t124 ^ 2;
t230 = t112 * t134;
t229 = t115 * t112;
t227 = t119 * t134;
t226 = t119 * t159;
t225 = t120 * t131;
t224 = t120 * t134;
t223 = t120 * t161;
t221 = t147 * t150;
t218 = t153 * t157;
t217 = t153 * t161;
t216 = t157 * t160;
t192 = t146 * t217;
t177 = t131 * t192 + t160;
t110 = t177 * t121;
t211 = -t110 + t160;
t129 = t134 ^ 2;
t109 = t112 * t129 + 0.1e1;
t202 = 0.2e1 / t109 ^ 2 * (-t129 * t232 - t134 * t229);
t201 = 0.2e1 * t232;
t200 = -0.2e1 * t231;
t116 = (-qJD(4) * t161 + qJD(1)) * t212 + (-t186 + (-qJD(1) * t161 + qJD(4)) * t160) * t150;
t135 = t160 * t149 + t161 * t219;
t130 = t135 ^ 2;
t127 = t130 * t218 + 0.1e1;
t158 = t156 * t157;
t199 = 0.2e1 * (t116 * t135 * t218 + (t153 * t158 * t210 - t157 * t188) * t130) / t127 ^ 2;
t198 = t152 * t231;
t197 = t112 * t227;
t194 = t131 * t222;
t191 = t156 * t217;
t190 = t157 * t210;
t182 = t111 * t202;
t181 = t112 * t202;
t180 = t152 * t199;
t178 = t146 * t198;
t176 = t134 * t201 + t229;
t133 = t150 * t213 - t212;
t175 = t131 * t221 - t133 * t146;
t174 = t133 * t156 - t135 * t216;
t125 = 0.1e1 / t127;
t118 = t135 * qJD(1) - t150 * t187 - t161 * t184 - t183;
t107 = 0.1e1 / t109;
t106 = t175 * t152 * t121;
t102 = (-t119 + (t120 * t194 + t119) * t121) * t134;
t101 = -t110 * t225 + (t211 * t226 + t223) * t149;
t100 = t120 * t150 * t159 - t119 * t133 + (-t119 * t215 - t225) * t106;
t98 = t177 * t200 + (t117 * t192 + t209 + (-t147 * t161 * t185 - t146 * t234) * t131) * t121;
t96 = -0.2e1 * t175 * t198 + (-t175 * t189 + (t117 * t221 - t118 * t146 + (t133 * t221 + (-0.2e1 * t148 * t150 ^ 2 - t146) * t131) * qJD(4)) * t152) * t121;
t1 = [t233 * t121 + 0.2e1 * t134 * t178, t98, 0, t96, 0, 0; t131 * t182 + (-t117 * t111 + (t102 * t115 + t131 * t99) * t112) * t107 + (t102 * t181 + (t102 * t201 + (t115 * t121 - t115 - (-t103 * t121 * t194 + t200) * t134) * t112 * t119 + (-(-0.2e1 * t131 * t178 - t103) * t230 + (-(t103 + t195) * t134 + t233 * t131) * t112 * t121) * t120) * t107) * t134, t101 * t134 * t181 + (-(-t98 * t225 + (t103 * t228 - t117 * t120) * t110) * t230 + t176 * t101 + (-t111 * t235 - (-t110 * t226 + t119 * t214 + t223) * t230) * t206) * t107 + (t182 * t235 + ((-t111 * t207 - (t211 * qJD(2) - t103) * t197) * t161 + (t111 * t210 + (t162 * t99 - (-t98 + t209) * t227 - (t211 * t103 - qJD(2)) * t224) * t112) * t159) * t107) * t149, 0 (t100 * t230 - t111 * t135) * t202 + (t116 * t111 + t176 * t100 - (-t118 + (-t103 * t150 - t149 * t96) * t159 - t171 * t106) * t197 + (-t135 * t99 - (t150 * t208 - t149 * t205 - t106 * t117 - t131 * t96 + (-t106 * t215 - t133) * t103) * t224) * t112) * t107, 0, 0; t174 * t180 + (t174 * t189 + (t116 * t216 - t118 * t156 + (-t133 * t216 + (0.2e1 * t158 * t160 ^ 2 + t156) * t135) * qJD(1)) * t152) * t125 (t135 * t191 + t150) * t199 + (-t116 * t191 + qJD(4) * t149 + (t156 * t234 - t190 * t217) * t135) * t125, 0, t134 * t156 * t180 + (t115 * t152 * t156 + (-t152 * t190 + t156 * t189) * t134) * t125, 0, 0;];
JaD_rot  = t1;
