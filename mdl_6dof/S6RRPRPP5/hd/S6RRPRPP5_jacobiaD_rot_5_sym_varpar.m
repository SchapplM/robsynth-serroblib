% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_jacobiaD_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:06
% EndTime: 2019-02-26 21:37:08
% DurationCPUTime: 1.03s
% Computational Cost: add. (1824->118), mult. (6168->263), div. (1114->15), fcn. (7752->9), ass. (0->114)
t154 = cos(qJ(2));
t155 = cos(qJ(1));
t206 = t154 * t155;
t231 = 0.2e1 * t206;
t151 = sin(qJ(2));
t150 = sin(qJ(4));
t205 = t155 * t150;
t152 = sin(qJ(1));
t153 = cos(qJ(4));
t208 = t152 * t153;
t134 = t151 * t208 + t205;
t207 = t154 * t153;
t124 = atan2(t134, t207);
t118 = sin(t124);
t119 = cos(t124);
t131 = t134 ^ 2;
t143 = 0.1e1 / t153 ^ 2;
t147 = 0.1e1 / t154 ^ 2;
t212 = t143 * t147;
t126 = t131 * t212 + 0.1e1;
t122 = 0.1e1 / t126;
t142 = 0.1e1 / t153;
t185 = t142 * t147 * t151;
t168 = t134 * t185 + t152;
t109 = t168 * t122;
t203 = t109 - t152;
t230 = t203 * t154 * t118 + t119 * t151;
t133 = t151 * t205 + t208;
t129 = 0.1e1 / t133 ^ 2;
t145 = t154 ^ 2;
t149 = t155 ^ 2;
t210 = t145 * t149;
t184 = t129 * t210;
t125 = 0.1e1 + t184;
t173 = qJD(4) * t151 + qJD(1);
t169 = t173 * t153;
t172 = qJD(1) * t151 + qJD(4);
t198 = qJD(2) * t155;
t178 = t154 * t198;
t117 = t155 * t169 + (-t172 * t152 + t178) * t150;
t128 = 0.1e1 / t133;
t220 = t117 * t128 * t129;
t171 = t210 * t220;
t199 = qJD(2) * t154;
t180 = t149 * t199;
t201 = qJD(1) * t155;
t182 = t152 * t201;
t229 = (-t171 + (-t145 * t182 - t151 * t180) * t129) / t125 ^ 2;
t146 = 0.1e1 / t154;
t204 = t155 * t153;
t209 = t152 * t150;
t135 = -t151 * t209 + t204;
t211 = t143 * t150;
t166 = t134 * t211 + t135 * t142;
t228 = t146 * t166;
t218 = t118 * t134;
t113 = t119 * t207 + t218;
t110 = 0.1e1 / t113;
t111 = 0.1e1 / t113 ^ 2;
t226 = -0.2e1 * t150;
t183 = t151 * t204;
t132 = -t183 + t209;
t127 = t132 ^ 2;
t108 = t127 * t111 + 0.1e1;
t116 = t134 * qJD(1) + t133 * qJD(4) - t153 * t178;
t221 = t116 * t111;
t195 = qJD(4) * t154;
t200 = qJD(2) * t151;
t165 = -t150 * t195 - t153 * t200;
t186 = t134 * t212;
t179 = t152 * t199;
t196 = qJD(4) * t153;
t114 = -qJD(1) * t183 - t153 * t179 - t155 * t196 + t173 * t209;
t213 = t142 * t146;
t188 = t114 * t213;
t102 = (-t165 * t186 - t188) * t122;
t164 = -t102 * t134 - t165;
t98 = (-t102 * t207 - t114) * t118 - t164 * t119;
t224 = t110 * t111 * t98;
t225 = 0.1e1 / t108 ^ 2 * (-t127 * t224 + t132 * t221);
t144 = t142 * t143;
t148 = t146 / t145;
t197 = qJD(4) * t150;
t177 = t147 * t197;
t223 = (-t114 * t186 + (t143 * t148 * t200 + t144 * t177) * t131) / t126 ^ 2;
t222 = t111 * t132;
t219 = t118 * t132;
t216 = t119 * t132;
t215 = t119 * t134;
t202 = qJD(1) * t152;
t194 = 0.2e1 * t225;
t193 = 0.2e1 * t224;
t192 = 0.2e1 * t229;
t191 = -0.2e1 * t223;
t190 = t110 * t225;
t189 = t111 * t219;
t187 = t134 * t213;
t181 = t147 * t200;
t176 = t111 * t194;
t175 = t132 * t193;
t174 = t132 * t231;
t170 = 0.2e1 * t213 * t223;
t167 = -t129 * t135 * t155 - t128 * t152;
t163 = t116 * t213 - (-t143 * t146 * t197 - t142 * t181) * t132;
t120 = 0.1e1 / t125;
t115 = t152 * t169 + (t172 * t155 + t179) * t150;
t106 = 0.1e1 / t108;
t105 = t122 * t228;
t101 = (-t118 + (-t119 * t187 + t118) * t122) * t132;
t100 = t109 * t215 - t230 * t153;
t99 = -t119 * t154 * t150 + t118 * t135 + (-t118 * t207 + t215) * t105;
t97 = t168 * t191 + (-t114 * t185 + t201 + (t143 * t151 * t177 + (0.2e1 * t148 * t151 ^ 2 + t146) * t142 * qJD(2)) * t134) * t122;
t95 = t191 * t228 + (t166 * t181 + (-t114 * t211 - t115 * t142 + (t135 * t211 + (0.2e1 * t144 * t150 ^ 2 + t142) * t134) * qJD(4)) * t146) * t122;
t1 = [-t163 * t122 + t132 * t170, t97, 0, t95, 0, 0; -0.2e1 * t134 * t190 + (-t114 * t110 + (-t101 * t116 - t134 * t98) * t111) * t106 + (t101 * t176 + (t101 * t193 + (-t116 * t122 + t116 - (t102 * t122 * t187 + t191) * t132) * t111 * t118 + (-(t134 * t170 - t102) * t222 + (-(t102 + t188) * t132 + t163 * t134) * t111 * t122) * t119) * t106) * t132, t100 * t132 * t176 + (-(t97 * t215 + (-t102 * t218 - t114 * t119) * t109) * t222 + (t175 - t221) * t100 + (t110 * t206 - t230 * t222) * t197) * t106 + (t190 * t231 + ((t110 * t198 - (t203 * qJD(2) + t102) * t189) * t151 + (t110 * t202 + (t155 * t98 - (-t97 + t201) * t219 - (-t203 * t102 - qJD(2)) * t216) * t111) * t154) * t106) * t153, 0 (-t110 * t133 + t99 * t222) * t194 + (t99 * t175 + t117 * t110 - (-t115 + (t102 * t150 - t153 * t95) * t154 + t164 * t105) * t189 + (-t99 * t116 - t133 * t98 - (t150 * t200 - t153 * t195 - t105 * t114 + t134 * t95 + (-t105 * t207 + t135) * t102) * t216) * t111) * t106, 0, 0; t167 * t154 * t192 + (t167 * t200 + ((qJD(1) * t128 - 0.2e1 * t135 * t220) * t155 + (-t115 * t155 + (-qJD(1) * t135 - t117) * t152) * t129) * t154) * t120 (-t128 * t151 * t155 - t150 * t184) * t192 + (t171 * t226 + (-t151 * t202 + t178) * t128 + ((-t117 * t155 + t180 * t226) * t151 + (t149 * t196 + t182 * t226) * t145) * t129) * t120, 0, t129 * t174 * t229 + (t174 * t220 + (-t116 * t206 + (t151 * t198 + t154 * t202) * t132) * t129) * t120, 0, 0;];
JaD_rot  = t1;
