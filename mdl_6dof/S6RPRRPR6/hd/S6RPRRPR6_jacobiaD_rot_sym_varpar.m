% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR6
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
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RPRRPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:32
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:13
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (2270->94), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->93)
	t136 = sin(qJ(1));
	t133 = t136 ^ 2;
	t132 = pkin(10) + qJ(3);
	t130 = sin(t132);
	t126 = t130 ^ 2;
	t131 = cos(t132);
	t128 = 0.1e1 / t131 ^ 2;
	t185 = t126 * t128;
	t121 = t133 * t185 + 0.1e1;
	t125 = t130 * t126;
	t127 = 0.1e1 / t131;
	t182 = t127 * t130;
	t146 = qJD(3) * (t125 * t127 * t128 + t182);
	t138 = cos(qJ(1));
	t174 = qJD(1) * t138;
	t183 = t126 * t136;
	t151 = t174 * t183;
	t191 = (t128 * t151 + t133 * t146) / t121 ^ 2;
	t201 = -0.2e1 * t191;
	t158 = 0.1e1 + t185;
	t200 = t136 * t158;
	t137 = cos(qJ(4));
	t176 = t137 * t138;
	t135 = sin(qJ(4));
	t178 = t136 * t135;
	t117 = t131 * t176 + t178;
	t179 = t136 * t130;
	t118 = atan2(-t179, -t131);
	t113 = cos(t118);
	t112 = sin(t118);
	t164 = t112 * t179;
	t102 = -t113 * t131 - t164;
	t99 = 0.1e1 / t102;
	t109 = 0.1e1 / t117;
	t100 = 0.1e1 / t102 ^ 2;
	t110 = 0.1e1 / t117 ^ 2;
	t119 = 0.1e1 / t121;
	t199 = t119 - 0.1e1;
	t134 = t138 ^ 2;
	t173 = qJD(3) * t131;
	t184 = t126 * t134;
	t172 = qJD(3) * t136;
	t160 = t128 * t172;
	t161 = t130 * t174;
	t93 = (-(-t131 * t172 - t161) * t127 + t126 * t160) * t119;
	t155 = t93 - t172;
	t156 = -t136 * t93 + qJD(3);
	t187 = t113 * t130;
	t88 = t156 * t187 + (t155 * t131 - t161) * t112;
	t196 = t99 * t100 * t88;
	t96 = t100 * t184 + 0.1e1;
	t198 = (-t184 * t196 + (t130 * t134 * t173 - t151) * t100) / t96 ^ 2;
	t94 = 0.1e1 / t96;
	t197 = t100 * t94;
	t177 = t136 * t137;
	t180 = t135 * t138;
	t116 = t131 * t180 - t177;
	t108 = t116 ^ 2;
	t107 = t108 * t110 + 0.1e1;
	t189 = t110 * t116;
	t153 = -qJD(1) * t131 + qJD(4);
	t154 = qJD(4) * t131 - qJD(1);
	t171 = qJD(3) * t138;
	t159 = t130 * t171;
	t98 = -t154 * t180 + (t153 * t136 - t159) * t137;
	t194 = t109 * t110 * t98;
	t147 = t131 * t178 + t176;
	t97 = t147 * qJD(1) - t117 * qJD(4) + t135 * t159;
	t195 = 0.1e1 / t107 ^ 2 * (-t108 * t194 - t97 * t189);
	t190 = t109 * t135;
	t188 = t112 * t131;
	t186 = t116 * t137;
	t181 = t130 * t138;
	t175 = qJD(1) * t136;
	t170 = 0.2e1 * t198;
	t169 = 0.2e1 * t196;
	t168 = -0.2e1 * t195;
	t167 = t99 * t198;
	t166 = t116 * t194;
	t165 = t94 * t173;
	t163 = t119 * t126 * t127;
	t157 = t127 * t201;
	t152 = t136 * t163;
	t150 = t158 * t138;
	t149 = t153 * t138;
	t148 = t110 * t186 - t190;
	t115 = -t131 * t177 + t180;
	t105 = 0.1e1 / t107;
	t104 = t119 * t200;
	t92 = (t199 * t130 * t112 - t113 * t152) * t138;
	t90 = -t136 * t188 + t187 + (-t113 * t179 + t188) * t104;
	t89 = t200 * t201 + (qJD(1) * t150 + 0.2e1 * t136 * t146) * t119;
	t1 = [t157 * t181 + (qJD(3) * t150 - t175 * t182) * t119, 0, t89, 0, 0, 0; (-t99 * t165 + (0.2e1 * t167 + (qJD(1) * t92 + t88) * t197) * t130) * t136 + ((-t92 * t165 + (t92 * t170 + ((0.2e1 * t130 * t191 - t93 * t152 - t199 * t173) * t112 + (t157 * t183 + t130 * t93 + (t125 * t160 - (t93 - 0.2e1 * t172) * t130) * t119) * t113) * t94 * t138) * t130) * t100 + (t92 * t169 + (-t99 + ((-t133 + t134) * t113 * t163 + t199 * t164) * t100) * qJD(1)) * t130 * t94) * t138, 0, (-t99 * t94 * t175 + (-0.2e1 * t167 + (-qJD(3) * t90 - t88) * t197) * t138) * t131 + (((-qJD(3) * t99 + t90 * t169) * t138 + (t90 * t175 + (-(-t104 * t174 - t136 * t89) * t113 - ((t104 * t136 - 0.1e1) * t93 + (-t104 + t136) * qJD(3)) * t112) * t181) * t100) * t94 + (t90 * t170 - ((t89 - t174) * t112 + (t155 * t104 + t156) * t113) * t94 * t131) * t100 * t138) * t130, 0, 0, 0; 0.2e1 * (t109 * t147 + t115 * t189) * t195 + (0.2e1 * t115 * t166 - t154 * t109 * t177 + (t130 * t172 + t149) * t190 + (t147 * t98 + t115 * t97 - t149 * t186 - (qJD(3) * t130 * t137 + t154 * t135) * t116 * t136) * t110) * t105, 0, t148 * t168 * t181 + (t148 * t131 * t171 + (-t148 * t175 + ((-qJD(4) * t109 - 0.2e1 * t166) * t137 + (-t137 * t97 + (-qJD(4) * t116 + t98) * t135) * t110) * t138) * t130) * t105, t168 + 0.2e1 * (-t105 * t110 * t97 + (-t105 * t194 - t110 * t195) * t116) * t116, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:13
	% DurationCPUTime: 1.07s
	% Computational Cost: add. (2618->95), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->94)
	t145 = pkin(10) + qJ(3);
	t141 = sin(t145);
	t137 = t141 ^ 2;
	t143 = cos(t145);
	t139 = 0.1e1 / t143 ^ 2;
	t195 = t137 * t139;
	t149 = sin(qJ(1));
	t212 = 0.2e1 * t149;
	t211 = t141 * t195;
	t147 = t149 ^ 2;
	t132 = t147 * t195 + 0.1e1;
	t130 = 0.1e1 / t132;
	t138 = 0.1e1 / t143;
	t150 = cos(qJ(1));
	t184 = qJD(1) * t150;
	t172 = t141 * t184;
	t182 = qJD(3) * t149;
	t104 = (-(-t143 * t182 - t172) * t138 + t182 * t195) * t130;
	t210 = t104 - t182;
	t146 = qJ(4) + pkin(11);
	t142 = sin(t146);
	t187 = t149 * t142;
	t144 = cos(t146);
	t189 = t144 * t150;
	t126 = t143 * t189 + t187;
	t188 = t149 * t141;
	t129 = atan2(-t188, -t143);
	t128 = cos(t129);
	t127 = sin(t129);
	t175 = t127 * t188;
	t113 = -t128 * t143 - t175;
	t110 = 0.1e1 / t113;
	t120 = 0.1e1 / t126;
	t111 = 0.1e1 / t113 ^ 2;
	t121 = 0.1e1 / t126 ^ 2;
	t209 = -0.2e1 * t141;
	t208 = t130 - 0.1e1;
	t197 = t128 * t141;
	t99 = (-t104 * t149 + qJD(3)) * t197 + (t143 * t210 - t172) * t127;
	t207 = t110 * t111 * t99;
	t160 = t143 * t187 + t189;
	t181 = qJD(3) * t150;
	t171 = t141 * t181;
	t105 = t160 * qJD(1) - qJD(4) * t126 + t142 * t171;
	t186 = t149 * t144;
	t190 = t143 * t150;
	t125 = t142 * t190 - t186;
	t119 = t125 ^ 2;
	t118 = t119 * t121 + 0.1e1;
	t200 = t121 * t125;
	t165 = -qJD(1) * t143 + qJD(4);
	t166 = qJD(4) * t143 - qJD(1);
	t191 = t142 * t150;
	t106 = -t166 * t191 + (t165 * t149 - t171) * t144;
	t204 = t106 * t120 * t121;
	t206 = (-t105 * t200 - t119 * t204) / t118 ^ 2;
	t205 = t104 * t141;
	t203 = t111 * t141;
	t193 = t138 * t141;
	t159 = qJD(3) * (t138 * t211 + t193);
	t163 = t137 * t149 * t184;
	t202 = (t139 * t163 + t147 * t159) / t132 ^ 2;
	t201 = t120 * t142;
	t199 = t125 * t144;
	t198 = t127 * t149;
	t196 = t137 * t138;
	t148 = t150 ^ 2;
	t194 = t137 * t148;
	t192 = t141 * t150;
	t185 = qJD(1) * t149;
	t183 = qJD(3) * t143;
	t109 = t111 * t194 + 0.1e1;
	t180 = 0.2e1 / t109 ^ 2 * (-t194 * t207 + (t141 * t148 * t183 - t163) * t111);
	t179 = 0.2e1 * t207;
	t178 = -0.2e1 * t206;
	t177 = t125 * t204;
	t176 = t111 * t192;
	t174 = t130 * t196;
	t170 = 0.1e1 + t195;
	t169 = t141 * t180;
	t168 = t202 * t209;
	t167 = t202 * t212;
	t164 = t149 * t174;
	t162 = t170 * t150;
	t161 = t121 * t199 - t201;
	t158 = t141 * t182 + t165 * t150;
	t124 = -t143 * t186 + t191;
	t117 = t170 * t149 * t130;
	t115 = 0.1e1 / t118;
	t107 = 0.1e1 / t109;
	t103 = (t208 * t141 * t127 - t128 * t164) * t150;
	t102 = -t143 * t198 + t197 + (t127 * t143 - t128 * t188) * t117;
	t100 = -t170 * t167 + (qJD(1) * t162 + t159 * t212) * t130;
	t1 = [t138 * t150 * t168 + (qJD(3) * t162 - t185 * t193) * t130, 0, t100, 0, 0, 0; (t110 * t169 + (-t110 * t183 + (qJD(1) * t103 + t99) * t203) * t107) * t149 + (t111 * t169 * t103 + (-((t104 * t164 + t208 * t183 + t168) * t127 + (t167 * t196 - t205 + (t205 + (t209 - t211) * t182) * t130) * t128) * t176 + (-t111 * t183 + t141 * t179) * t103 + (-t110 + ((-t147 + t148) * t128 * t174 + t208 * t175) * t111) * t141 * qJD(1)) * t107) * t150, 0, (t102 * t203 - t110 * t143) * t150 * t180 + ((-t110 * t185 + (-qJD(3) * t102 - t99) * t150 * t111) * t143 + (-t110 * t181 - (-t100 * t128 * t149 - t210 * t127 + (-qJD(3) * t127 + t104 * t198 - t128 * t184) * t117) * t176 + (t111 * t185 + t150 * t179) * t102 - ((t100 - t184) * t127 + ((-t117 * t149 + 0.1e1) * qJD(3) + (t117 - t149) * t104) * t128) * t111 * t190) * t141) * t107, 0, 0, 0; 0.2e1 * (t120 * t160 + t124 * t200) * t206 + (0.2e1 * t124 * t177 - t166 * t120 * t186 + t158 * t201 + (-t166 * t125 * t187 + t124 * t105 + t106 * t160 - t158 * t199) * t121) * t115, 0, t161 * t178 * t192 + (t161 * t143 * t181 + (-t161 * t185 + ((-qJD(4) * t120 - 0.2e1 * t177) * t144 + (-t105 * t144 + (-qJD(4) * t125 + t106) * t142) * t121) * t150) * t141) * t115, t178 + 0.2e1 * (-t105 * t115 * t121 + (-t115 * t204 - t121 * t206) * t125) * t125, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:32:12
	% EndTime: 2019-10-10 01:32:13
	% DurationCPUTime: 1.10s
	% Computational Cost: add. (3326->96), mult. (2734->203), div. (498->12), fcn. (3199->9), ass. (0->95)
	t162 = pkin(10) + qJ(3);
	t159 = sin(t162);
	t155 = t159 ^ 2;
	t160 = cos(t162);
	t157 = 0.1e1 / t160 ^ 2;
	t210 = t155 * t157;
	t166 = sin(qJ(1));
	t228 = 0.2e1 * t166;
	t227 = t159 * t210;
	t164 = t166 ^ 2;
	t149 = t164 * t210 + 0.1e1;
	t147 = 0.1e1 / t149;
	t156 = 0.1e1 / t160;
	t167 = cos(qJ(1));
	t201 = qJD(1) * t167;
	t189 = t159 * t201;
	t199 = qJD(3) * t166;
	t122 = (-(-t160 * t199 - t189) * t156 + t199 * t210) * t147;
	t226 = t122 - t199;
	t161 = qJ(4) + pkin(11) + qJ(6);
	t153 = cos(t161);
	t152 = sin(t161);
	t205 = t166 * t152;
	t206 = t160 * t167;
	t142 = t153 * t206 + t205;
	t203 = t166 * t159;
	t146 = atan2(-t203, -t160);
	t145 = cos(t146);
	t144 = sin(t146);
	t192 = t144 * t203;
	t132 = -t145 * t160 - t192;
	t129 = 0.1e1 / t132;
	t136 = 0.1e1 / t142;
	t130 = 0.1e1 / t132 ^ 2;
	t137 = 0.1e1 / t142 ^ 2;
	t225 = -0.2e1 * t159;
	t224 = t147 - 0.1e1;
	t213 = t145 * t159;
	t115 = (-t122 * t166 + qJD(3)) * t213 + (t160 * t226 - t189) * t144;
	t223 = t115 * t129 * t130;
	t163 = qJD(4) + qJD(6);
	t177 = t153 * t167 + t160 * t205;
	t198 = qJD(3) * t167;
	t188 = t159 * t198;
	t120 = t177 * qJD(1) - t142 * t163 + t152 * t188;
	t204 = t166 * t153;
	t141 = t152 * t206 - t204;
	t135 = t141 ^ 2;
	t128 = t135 * t137 + 0.1e1;
	t216 = t137 * t141;
	t182 = -qJD(1) * t160 + t163;
	t183 = t160 * t163 - qJD(1);
	t212 = t152 * t167;
	t121 = -t183 * t212 + (t182 * t166 - t188) * t153;
	t221 = t121 * t136 * t137;
	t222 = (-t120 * t216 - t135 * t221) / t128 ^ 2;
	t220 = t122 * t159;
	t219 = t130 * t159;
	t208 = t156 * t159;
	t176 = qJD(3) * (t156 * t227 + t208);
	t180 = t155 * t166 * t201;
	t218 = (t157 * t180 + t164 * t176) / t149 ^ 2;
	t217 = t136 * t152;
	t215 = t141 * t153;
	t214 = t144 * t166;
	t211 = t155 * t156;
	t165 = t167 ^ 2;
	t209 = t155 * t165;
	t207 = t159 * t167;
	t202 = qJD(1) * t166;
	t200 = qJD(3) * t160;
	t125 = t130 * t209 + 0.1e1;
	t197 = 0.2e1 * (-t209 * t223 + (t159 * t165 * t200 - t180) * t130) / t125 ^ 2;
	t196 = 0.2e1 * t223;
	t195 = -0.2e1 * t222;
	t194 = t141 * t221;
	t193 = t130 * t207;
	t191 = t147 * t211;
	t187 = 0.1e1 + t210;
	t186 = t159 * t197;
	t185 = t218 * t225;
	t184 = t218 * t228;
	t181 = t166 * t191;
	t179 = t187 * t167;
	t178 = t137 * t215 - t217;
	t175 = t159 * t199 + t182 * t167;
	t140 = -t160 * t204 + t212;
	t134 = t187 * t166 * t147;
	t126 = 0.1e1 / t128;
	t123 = 0.1e1 / t125;
	t119 = (t224 * t159 * t144 - t145 * t181) * t167;
	t118 = -t160 * t214 + t213 + (t144 * t160 - t145 * t203) * t134;
	t116 = -t187 * t184 + (qJD(1) * t179 + t176 * t228) * t147;
	t113 = t195 + 0.2e1 * (-t120 * t126 * t137 + (-t126 * t221 - t137 * t222) * t141) * t141;
	t1 = [t156 * t167 * t185 + (qJD(3) * t179 - t202 * t208) * t147, 0, t116, 0, 0, 0; (t129 * t186 + (-t129 * t200 + (qJD(1) * t119 + t115) * t219) * t123) * t166 + (t130 * t186 * t119 + (-((t122 * t181 + t224 * t200 + t185) * t144 + (t184 * t211 - t220 + (t220 + (t225 - t227) * t199) * t147) * t145) * t193 + (-t130 * t200 + t159 * t196) * t119 + (-t129 + ((-t164 + t165) * t145 * t191 + t224 * t192) * t130) * t159 * qJD(1)) * t123) * t167, 0, (t118 * t219 - t129 * t160) * t167 * t197 + ((-t129 * t202 + (-qJD(3) * t118 - t115) * t167 * t130) * t160 + (-t129 * t198 - (-t116 * t145 * t166 - t226 * t144 + (-qJD(3) * t144 + t122 * t214 - t145 * t201) * t134) * t193 + (t130 * t202 + t167 * t196) * t118 - ((t116 - t201) * t144 + ((-t134 * t166 + 0.1e1) * qJD(3) + (t134 - t166) * t122) * t145) * t130 * t206) * t159) * t123, 0, 0, 0; 0.2e1 * (t136 * t177 + t140 * t216) * t222 + (0.2e1 * t140 * t194 - t183 * t136 * t204 + t175 * t217 + (-t183 * t141 * t205 + t140 * t120 + t121 * t177 - t175 * t215) * t137) * t126, 0, t178 * t195 * t207 + (t178 * t160 * t198 + (-t178 * t202 + ((-t136 * t163 - 0.2e1 * t194) * t153 + (-t120 * t153 + (-t141 * t163 + t121) * t152) * t137) * t167) * t159) * t126, t113, 0, t113;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end