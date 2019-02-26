% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_jacobiaD_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:33:08
% EndTime: 2019-02-26 20:33:09
% DurationCPUTime: 1.14s
% Computational Cost: add. (1824->118), mult. (6168->266), div. (1114->15), fcn. (7752->9), ass. (0->111)
t153 = sin(qJ(4));
t152 = sin(qJ(5));
t154 = sin(qJ(1));
t210 = t154 * t152;
t184 = t153 * t210;
t155 = cos(qJ(5));
t157 = cos(qJ(1));
t205 = t157 * t155;
t133 = t184 - t205;
t156 = cos(qJ(4));
t208 = t156 * t152;
t125 = atan2(-t133, t208);
t119 = sin(t125);
t120 = cos(t125);
t128 = t133 ^ 2;
t145 = 0.1e1 / t152 ^ 2;
t149 = 0.1e1 / t156 ^ 2;
t213 = t145 * t149;
t127 = t128 * t213 + 0.1e1;
t123 = 0.1e1 / t127;
t144 = 0.1e1 / t152;
t186 = t144 * t149 * t153;
t172 = t133 * t186 + t154;
t110 = t172 * t123;
t230 = t110 - t154;
t232 = t230 * t119 * t156 - t120 * t153;
t207 = t156 * t157;
t231 = -0.2e1 * t207;
t175 = qJD(1) * t153 + qJD(5);
t198 = qJD(5) * t155;
t179 = t153 * t198;
t199 = qJD(4) * t157;
t180 = t156 * t199;
t202 = qJD(1) * t157;
t183 = t155 * t202;
t115 = -t152 * t180 - t157 * t179 + t175 * t210 - t183;
t206 = t157 * t152;
t209 = t154 * t155;
t136 = t153 * t206 + t209;
t148 = 0.1e1 / t156;
t201 = qJD(4) * t153;
t182 = t149 * t201;
t214 = t144 * t148;
t229 = -(-t145 * t148 * t198 + t144 * t182) * t136 + t115 * t214;
t220 = t119 * t133;
t114 = t120 * t208 - t220;
t111 = 0.1e1 / t114;
t137 = t153 * t205 - t210;
t130 = 0.1e1 / t137;
t112 = 0.1e1 / t114 ^ 2;
t131 = 0.1e1 / t137 ^ 2;
t228 = 0.2e1 * t155;
t129 = t136 ^ 2;
t109 = t129 * t112 + 0.1e1;
t222 = t115 * t112;
t197 = qJD(5) * t156;
t168 = -t152 * t201 + t155 * t197;
t187 = t133 * t213;
t200 = qJD(4) * t156;
t166 = t154 * t200 + t175 * t157;
t176 = qJD(5) * t153 + qJD(1);
t117 = t166 * t152 + t176 * t209;
t189 = t117 * t214;
t103 = (t168 * t187 - t189) * t123;
t165 = t103 * t133 - t168;
t99 = (-t103 * t208 - t117) * t119 - t165 * t120;
t226 = t111 * t112 * t99;
t227 = 0.1e1 / t109 ^ 2 * (-t129 * t226 - t136 * t222);
t147 = t156 ^ 2;
t151 = t157 ^ 2;
t211 = t147 * t151;
t185 = t131 * t211;
t126 = 0.1e1 + t185;
t116 = -t176 * t206 + (-t175 * t154 + t180) * t155;
t221 = t116 * t130 * t131;
t173 = t211 * t221;
t181 = t151 * t200;
t225 = (-t173 + (-t147 * t154 * t202 - t153 * t181) * t131) / t126 ^ 2;
t146 = t144 * t145;
t150 = t148 / t147;
t224 = (t117 * t187 + (t145 * t150 * t201 - t146 * t149 * t198) * t128) / t127 ^ 2;
t223 = t112 * t136;
t219 = t119 * t136;
t217 = t120 * t133;
t216 = t120 * t136;
t212 = t145 * t155;
t203 = qJD(1) * t154;
t196 = 0.2e1 * t227;
t195 = 0.2e1 * t226;
t194 = 0.2e1 * t225;
t193 = t111 * t227;
t192 = t148 * t224;
t191 = t112 * t219;
t188 = t133 * t214;
t178 = t112 * t196;
t177 = t136 * t231;
t174 = t144 * t192;
t171 = t136 * t195 + t222;
t135 = t153 * t209 + t206;
t170 = -t131 * t135 * t157 + t130 * t154;
t169 = t133 * t212 - t135 * t144;
t121 = 0.1e1 / t126;
t118 = -qJD(5) * t184 - t152 * t203 + t166 * t155;
t107 = 0.1e1 / t109;
t106 = t169 * t148 * t123;
t102 = (-t119 + (t120 * t188 + t119) * t123) * t136;
t101 = t110 * t217 + t232 * t152;
t100 = t120 * t156 * t155 - t119 * t135 + (-t119 * t208 - t217) * t106;
t98 = 0.2e1 * t172 * t224 + (-t117 * t186 - t202 + (t179 * t213 + (-0.2e1 * t150 * t153 ^ 2 - t148) * t144 * qJD(4)) * t133) * t123;
t96 = -0.2e1 * t169 * t192 + (t169 * t182 + (t117 * t212 - t118 * t144 + (t135 * t212 + (-0.2e1 * t146 * t155 ^ 2 - t144) * t133) * qJD(5)) * t148) * t123;
t1 = [t229 * t123 + 0.2e1 * t136 * t174, 0, 0, t98, t96, 0; 0.2e1 * t133 * t193 + (-t117 * t111 + (t102 * t115 + t133 * t99) * t112) * t107 + (t102 * t178 + (t102 * t195 + (t115 * t123 - t115 - (-t103 * t123 * t188 - 0.2e1 * t224) * t136) * t112 * t119 + (-(-0.2e1 * t133 * t174 - t103) * t223 + (-(t103 + t189) * t136 + t229 * t133) * t112 * t123) * t120) * t107) * t136, 0, 0, t101 * t136 * t178 + (-(-t98 * t217 - (t103 * t220 - t117 * t120) * t110) * t223 + t171 * t101 + (t111 * t207 - t232 * t223) * t198) * t107 + (t193 * t231 + ((-t111 * t199 - (-qJD(4) * t230 + t103) * t191) * t153 + (-t111 * t203 + (-t157 * t99 - (-t98 - t202) * t219 - (t103 * t230 - qJD(4)) * t216) * t112) * t156) * t107) * t152 (t100 * t223 - t111 * t137) * t196 + (t116 * t111 + t171 * t100 - (-t118 + (-t103 * t155 - t152 * t96) * t156 + t165 * t106) * t191 + (-t137 * t99 - (-t155 * t201 - t152 * t197 - t106 * t117 - t133 * t96 + (-t106 * t208 - t135) * t103) * t216) * t112) * t107, 0; t170 * t156 * t194 + (t170 * t201 + ((-qJD(1) * t130 - 0.2e1 * t135 * t221) * t157 + (t118 * t157 + (-qJD(1) * t135 + t116) * t154) * t131) * t156) * t121, 0, 0 (t130 * t153 * t157 + t155 * t185) * t194 + (t173 * t228 + (t153 * t203 - t180) * t130 + ((t116 * t157 + t181 * t228) * t153 + (qJD(5) * t151 * t152 + 0.2e1 * t154 * t183) * t147) * t131) * t121, t131 * t177 * t225 + (t177 * t221 + (-t115 * t207 + (-t153 * t199 - t156 * t203) * t136) * t131) * t121, 0;];
JaD_rot  = t1;
