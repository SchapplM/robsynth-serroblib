% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP8_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_jacobiaD_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:09
% EndTime: 2019-02-26 22:13:10
% DurationCPUTime: 1.00s
% Computational Cost: add. (1824->123), mult. (6168->270), div. (1114->15), fcn. (7752->9), ass. (0->112)
t155 = sin(qJ(3));
t158 = cos(qJ(2));
t157 = cos(qJ(3));
t159 = cos(qJ(1));
t211 = t159 * t157;
t230 = sin(qJ(1));
t136 = t230 * t155 + t158 * t211;
t130 = 0.1e1 / t136 ^ 2;
t156 = sin(qJ(2));
t150 = t156 ^ 2;
t154 = t159 ^ 2;
t215 = t150 * t154;
t196 = t130 * t215;
t125 = 0.1e1 + t196;
t184 = qJD(1) * t230;
t208 = qJD(2) * t159;
t188 = t156 * t208;
t169 = t158 * t184 + t188;
t183 = t230 * qJD(3);
t212 = t159 * t155;
t115 = (-qJD(3) * t158 + qJD(1)) * t212 + (t183 - t169) * t157;
t129 = 0.1e1 / t136;
t225 = t115 * t129 * t130;
t177 = t215 * t225;
t189 = qJD(2) * t154 * t156;
t233 = (-t177 + (-t150 * t159 * t184 + t158 * t189) * t130) / t125 ^ 2;
t213 = t156 * t159;
t191 = t230 * t158;
t132 = t155 * t191 + t211;
t174 = t155 * t183;
t205 = qJD(3) * t159;
t186 = t157 * t205;
t114 = t132 * qJD(1) + t155 * t188 - t158 * t186 - t174;
t135 = -t230 * t157 + t158 * t212;
t147 = 0.1e1 / t155;
t148 = 0.1e1 / t155 ^ 2;
t151 = 0.1e1 / t156;
t152 = 0.1e1 / t156 ^ 2;
t209 = qJD(2) * t158;
t190 = t152 * t209;
t206 = qJD(3) * t157;
t218 = t147 * t151;
t232 = (t148 * t151 * t206 + t147 * t190) * t135 + t114 * t218;
t214 = t156 * t155;
t124 = atan2(-t132, t214);
t119 = cos(t124);
t118 = sin(t124);
t224 = t118 * t132;
t113 = t119 * t214 - t224;
t110 = 0.1e1 / t113;
t111 = 0.1e1 / t113 ^ 2;
t231 = 0.2e1 * t135;
t127 = t132 ^ 2;
t217 = t148 * t152;
t126 = t127 * t217 + 0.1e1;
t122 = 0.1e1 / t126;
t170 = t155 * t209 + t156 * t206;
t194 = t132 * t217;
t192 = t230 * t156;
t175 = qJD(2) * t192;
t176 = t157 * t184;
t210 = qJD(1) * t159;
t116 = t157 * t183 * t158 - t176 + (t210 * t158 - t175 - t205) * t155;
t197 = t116 * t218;
t102 = (t170 * t194 - t197) * t122;
t167 = -t102 * t132 + t170;
t98 = (-t102 * t214 - t116) * t118 + t167 * t119;
t229 = t110 * t111 * t98;
t149 = t147 * t148;
t153 = t151 / t150;
t187 = t152 * t206;
t228 = (t116 * t194 + (-t148 * t153 * t209 - t149 * t187) * t127) / t126 ^ 2;
t227 = t111 * t135;
t226 = t114 * t111;
t223 = t118 * t135;
t222 = t118 * t156;
t221 = t119 * t132;
t220 = t119 * t135;
t219 = t119 * t158;
t216 = t148 * t157;
t207 = qJD(3) * t155;
t128 = t135 ^ 2;
t108 = t111 * t128 + 0.1e1;
t204 = 0.2e1 / t108 ^ 2 * (-t128 * t229 - t135 * t226);
t203 = 0.2e1 * t229;
t202 = 0.2e1 * t233;
t201 = -0.2e1 * t228;
t200 = t151 * t228;
t199 = t111 * t223;
t195 = t132 * t218;
t193 = t147 * t152 * t158;
t172 = t132 * t193 + t230;
t109 = t172 * t122;
t185 = t230 - t109;
t182 = t110 * t204;
t181 = t111 * t204;
t180 = t135 * t203;
t179 = t213 * t231;
t178 = t147 * t200;
t134 = t157 * t191 - t212;
t173 = t132 * t216 - t134 * t147;
t171 = t130 * t134 * t159 - t230 * t129;
t120 = 0.1e1 / t125;
t117 = t136 * qJD(1) - t157 * t175 - t158 * t174 - t186;
t106 = 0.1e1 / t108;
t105 = t173 * t151 * t122;
t101 = (-t118 + (t119 * t195 + t118) * t122) * t135;
t100 = -t109 * t221 + (t185 * t222 + t219) * t155;
t99 = t119 * t156 * t157 - t118 * t134 + (-t118 * t214 - t221) * t105;
t97 = t172 * t201 + (t116 * t193 + t210 + (-t148 * t158 * t187 + (-0.2e1 * t153 * t158 ^ 2 - t151) * t147 * qJD(2)) * t132) * t122;
t95 = -0.2e1 * t173 * t200 + (-t173 * t190 + (t116 * t216 - t117 * t147 + (t134 * t216 + (-0.2e1 * t149 * t157 ^ 2 - t147) * t132) * qJD(3)) * t151) * t122;
t1 = [t232 * t122 + t178 * t231, t97, t95, 0, 0, 0; t132 * t182 + (-t116 * t110 + (t101 * t114 + t132 * t98) * t111) * t106 + (t101 * t181 + (t101 * t203 + (t114 * t122 - t114 - (-t102 * t122 * t195 + t201) * t135) * t111 * t118 + (-(-0.2e1 * t132 * t178 - t102) * t227 + (-(t102 + t197) * t135 + t232 * t132) * t111 * t122) * t119) * t106) * t135, t100 * t135 * t181 + (-(-t97 * t221 + (t102 * t224 - t116 * t119) * t109) * t227 + (t180 + t226) * t100 + (-t110 * t213 - (-t109 * t222 + t118 * t192 + t219) * t227) * t206) * t106 + (t182 * t213 + ((-t110 * t208 - (t185 * qJD(2) - t102) * t199) * t158 + (t110 * t184 + (t159 * t98 - (-t97 + t210) * t223 - (t185 * t102 - qJD(2)) * t220) * t111) * t156) * t106) * t155 (-t110 * t136 + t99 * t227) * t204 + (t99 * t180 + t115 * t110 - (-t117 + (-t102 * t157 - t155 * t95) * t156 - t167 * t105) * t199 + (t99 * t114 - t136 * t98 - (t157 * t209 - t156 * t207 - t105 * t116 - t132 * t95 + (-t105 * t214 - t134) * t102) * t220) * t111) * t106, 0, 0, 0; t171 * t156 * t202 + (-t171 * t209 + ((qJD(1) * t129 + 0.2e1 * t134 * t225) * t159 + (-t230 * t115 - t117 * t159 + t134 * t184) * t130) * t156) * t120 (t129 * t158 * t159 + t157 * t196) * t202 + (0.2e1 * t157 * t177 + t169 * t129 + ((t115 * t159 - 0.2e1 * t157 * t189) * t158 + (t154 * t207 + 0.2e1 * t159 * t176) * t150) * t130) * t120, t130 * t179 * t233 + (t179 * t225 + (t114 * t213 + (t156 * t184 - t158 * t208) * t135) * t130) * t120, 0, 0, 0;];
JaD_rot  = t1;
