% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_jacobiaD_rot_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:59:51
% EndTime: 2019-02-26 20:59:52
% DurationCPUTime: 0.89s
% Computational Cost: add. (1824->120), mult. (6168->265), div. (1114->15), fcn. (7752->9), ass. (0->113)
t153 = sin(qJ(1));
t155 = cos(qJ(3));
t208 = t153 * t155;
t152 = sin(qJ(3));
t151 = sin(qJ(4));
t228 = cos(qJ(1));
t190 = t228 * t151;
t154 = cos(qJ(4));
t209 = t153 * t154;
t135 = t152 * t190 + t209;
t207 = t155 * t151;
t125 = atan2(t135, t207);
t120 = cos(t125);
t119 = sin(t125);
t219 = t119 * t135;
t114 = t120 * t207 + t219;
t111 = 0.1e1 / t114;
t134 = t152 * t209 + t190;
t129 = 0.1e1 / t134;
t143 = 0.1e1 / t151;
t148 = 0.1e1 / t155;
t112 = 0.1e1 / t114 ^ 2;
t130 = 0.1e1 / t134 ^ 2;
t144 = 0.1e1 / t151 ^ 2;
t189 = t228 * t154;
t210 = t153 * t151;
t133 = t152 * t210 - t189;
t128 = t133 ^ 2;
t109 = t112 * t128 + 0.1e1;
t182 = qJD(1) * t228;
t172 = t152 * t182;
t165 = t228 * qJD(4) + t172;
t177 = qJD(4) * t152 + qJD(1);
t203 = qJD(3) * t155;
t186 = t153 * t203;
t117 = t177 * t209 + (t165 + t186) * t151;
t222 = t117 * t112;
t132 = t135 ^ 2;
t149 = 0.1e1 / t155 ^ 2;
t213 = t144 * t149;
t127 = t132 * t213 + 0.1e1;
t123 = 0.1e1 / t127;
t201 = qJD(4) * t155;
t205 = qJD(3) * t152;
t166 = -t151 * t205 + t154 * t201;
t192 = t135 * t213;
t188 = t228 * t155;
t173 = qJD(3) * t188;
t174 = t154 * t182;
t175 = t152 * t189;
t115 = -qJD(4) * t175 - t151 * t173 - t174 + (qJD(1) * t152 + qJD(4)) * t210;
t214 = t143 * t148;
t195 = t115 * t214;
t103 = (-t166 * t192 - t195) * t123;
t164 = -t103 * t135 - t166;
t99 = (-t103 * t207 - t115) * t119 - t164 * t120;
t226 = t111 * t112 * t99;
t227 = 0.1e1 / t109 ^ 2 * (-t128 * t226 + t133 * t222);
t146 = t153 ^ 2;
t147 = t155 ^ 2;
t211 = t146 * t147;
t194 = t130 * t211;
t126 = 0.1e1 + t194;
t185 = t154 * t203;
t118 = t165 * t154 + (-t177 * t151 + t185) * t153;
t221 = t118 * t129 * t130;
t176 = t211 * t221;
t225 = (-t176 + (-t146 * t152 * t203 + t147 * t153 * t182) * t130) / t126 ^ 2;
t145 = t143 * t144;
t150 = t148 / t147;
t202 = qJD(4) * t154;
t184 = t149 * t202;
t224 = (-t115 * t192 + (t144 * t150 * t205 - t145 * t184) * t132) / t127 ^ 2;
t223 = t112 * t133;
t220 = t119 * t133;
t218 = t119 * t155;
t217 = t120 * t133;
t216 = t120 * t135;
t215 = t120 * t152;
t212 = t144 * t154;
t206 = qJD(1) * t153;
t204 = qJD(3) * t153;
t200 = 0.2e1 * t227;
t199 = 0.2e1 * t226;
t198 = 0.2e1 * t225;
t197 = -0.2e1 * t224;
t196 = t112 * t220;
t193 = t135 * t214;
t191 = t143 * t149 * t152;
t187 = t149 * t205;
t168 = t135 * t191 + t228;
t110 = t168 * t123;
t183 = t228 - t110;
t181 = -0.2e1 * t111 * t227;
t180 = t112 * t200;
t179 = 0.2e1 * t148 * t224;
t178 = -0.2e1 * t133 * t208;
t171 = t143 * t179;
t170 = t133 * t199 - t222;
t136 = t175 - t210;
t169 = t135 * t212 - t136 * t143;
t167 = t130 * t136 * t153 - t228 * t129;
t163 = t117 * t214 - (t144 * t148 * t202 - t143 * t187) * t133;
t121 = 0.1e1 / t126;
t116 = t134 * qJD(1) + qJD(4) * t135 - t154 * t173;
t107 = 0.1e1 / t109;
t106 = t169 * t148 * t123;
t102 = (-t119 + (-t120 * t193 + t119) * t123) * t133;
t101 = t110 * t216 + (t183 * t218 - t215) * t151;
t100 = t120 * t154 * t155 + t119 * t136 - (-t119 * t207 + t216) * t106;
t98 = t168 * t197 + (-t115 * t191 - t206 + (-t144 * t152 * t184 + (0.2e1 * t150 * t152 ^ 2 + t148) * t143 * qJD(3)) * t135) * t123;
t96 = t169 * t179 + (-t169 * t187 + (t115 * t212 - t116 * t143 + (-t136 * t212 + (0.2e1 * t145 * t154 ^ 2 + t143) * t135) * qJD(4)) * t148) * t123;
t1 = [-t163 * t123 + t133 * t171, 0, t98, t96, 0, 0; t135 * t181 + (-t115 * t111 + (-t102 * t117 - t135 * t99) * t112) * t107 + (t102 * t180 + (t102 * t199 + (-t117 * t123 + t117 - (t103 * t123 * t193 + t197) * t133) * t112 * t119 + (-(t135 * t171 - t103) * t223 + (-(t103 + t195) * t133 + t163 * t135) * t112 * t123) * t120) * t107) * t133, 0, t101 * t133 * t180 + (-(t98 * t216 + (-t103 * t219 - t115 * t120) * t110) * t223 + t170 * t101 + (t111 * t208 - (-t110 * t218 + t119 * t188 - t215) * t223) * t202) * t107 + (t181 * t208 + ((-t111 * t204 - (-t183 * qJD(3) + t103) * t196) * t152 + (t111 * t182 + (-t153 * t99 - (-t98 - t206) * t220 - (t183 * t103 - qJD(3)) * t217) * t112) * t155) * t107) * t151 (t100 * t223 - t111 * t134) * t200 + (t118 * t111 + t170 * t100 - (-t116 + (-t103 * t154 - t151 * t96) * t155 - t164 * t106) * t196 + (-t134 * t99 - (-t154 * t205 - t151 * t201 + t106 * t115 + t135 * t96 + (t106 * t207 + t136) * t103) * t217) * t112) * t107, 0, 0; t167 * t155 * t198 + (t167 * t205 + ((-qJD(1) * t129 + 0.2e1 * t136 * t221) * t153 + (t116 * t153 - t228 * t118 - t136 * t182) * t130) * t155) * t121, 0 (t129 * t152 * t153 + t154 * t194) * t198 + (0.2e1 * t154 * t176 + (-t172 - t186) * t129 + ((t118 * t152 - 0.2e1 * t147 * t174) * t153 + (qJD(4) * t147 * t151 + 0.2e1 * t152 * t185) * t146) * t130) * t121, t130 * t178 * t225 + (t178 * t221 + (t117 * t208 + (-t152 * t204 + t155 * t182) * t133) * t130) * t121, 0, 0;];
JaD_rot  = t1;
