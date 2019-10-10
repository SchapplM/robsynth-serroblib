% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR12
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
%   Wie in S6RRPRPR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR12_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR12_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (215->39), mult. (853->106), div. (126->12), fcn. (1047->9), ass. (0->54)
	t99 = sin(pkin(6));
	t93 = t99 ^ 2;
	t100 = cos(pkin(6));
	t95 = 0.1e1 / t100 ^ 2;
	t104 = cos(qJ(1));
	t98 = t104 ^ 2;
	t89 = t98 * t93 * t95 + 0.1e1;
	t102 = sin(qJ(1));
	t97 = t102 ^ 2;
	t126 = 0.1e1 / t89 ^ 2 * t97;
	t131 = t126 * t95;
	t122 = t104 * t99;
	t88 = atan2(t122, t100);
	t84 = sin(t88);
	t85 = cos(t88);
	t72 = t85 * t100 + t84 * t122;
	t67 = 0.1e1 / t72;
	t103 = cos(qJ(2));
	t118 = t104 * t103;
	t101 = sin(qJ(2));
	t121 = t102 * t101;
	t113 = t100 * t121 - t118;
	t77 = 0.1e1 / t113;
	t94 = 0.1e1 / t100;
	t68 = 0.1e1 / t72 ^ 2;
	t78 = 0.1e1 / t113 ^ 2;
	t119 = t104 * t101;
	t120 = t102 * t103;
	t81 = -t100 * t119 - t120;
	t82 = t100 * t120 + t119;
	t71 = t81 * qJD(1) - t82 * qJD(2);
	t128 = t71 * t77 * t78;
	t115 = t100 * t118;
	t70 = -qJD(1) * t115 - qJD(2) * t118 + (qJD(2) * t100 + qJD(1)) * t121;
	t129 = t70 * t78;
	t76 = t82 ^ 2;
	t75 = t76 * t78 + 0.1e1;
	t130 = (t76 * t128 - t82 * t129) / t75 ^ 2;
	t127 = t81 * t82;
	t125 = t93 * t94;
	t124 = t102 * t68;
	t123 = t104 * t68;
	t117 = qJD(1) * t104;
	t86 = 0.1e1 / t89;
	t116 = (t86 - 0.1e1) * t99;
	t114 = -0.2e1 * t94 * t131;
	t80 = t115 - t121;
	t63 = (-t104 * t85 * t86 * t125 + t84 * t116) * t102;
	t92 = t99 * t93;
	t73 = 0.1e1 / t75;
	t69 = t67 * t68;
	t66 = t97 * t93 * t68 + 0.1e1;
	t62 = qJD(1) * t63;
	t1 = [(-t86 * t94 * t99 + t92 * t114) * t117, 0, 0, 0, 0, 0; (0.2e1 * (-t104 * t67 + t63 * t124) / t66 ^ 2 * (-t62 * t69 * t97 + t117 * t124) * t93 + ((0.2e1 * t102 * t63 * t69 - t123) * t62 + (-t63 * t123 + (-t67 + (-t92 * t131 - t116) * t84 * t123 - (t93 ^ 2 * t98 * t114 + (-t126 + (0.2e1 * t97 - t98) * t86) * t125) * t68 * t85) * t102) * qJD(1)) / t66) * t99, 0, 0, 0, 0, 0; 0.2e1 * (t78 * t127 + t77 * t80) * t130 + (-(-t82 * qJD(1) + t81 * qJD(2)) * t77 - 0.2e1 * t127 * t128 + (-t80 * t71 - (t113 * qJD(1) - t80 * qJD(2)) * t82 + t81 * t70) * t78) * t73, -0.2e1 * t130 + 0.2e1 * (-t73 * t129 + (t73 * t128 - t78 * t130) * t82) * t82, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:32
	% DurationCPUTime: 0.77s
	% Computational Cost: add. (1114->72), mult. (3196->173), div. (656->14), fcn. (4222->9), ass. (0->74)
	t125 = cos(qJ(2));
	t126 = cos(qJ(1));
	t156 = cos(pkin(6));
	t141 = t126 * t156;
	t139 = t125 * t141;
	t123 = sin(qJ(2));
	t124 = sin(qJ(1));
	t152 = t124 * t123;
	t104 = -t139 + t152;
	t122 = sin(pkin(6));
	t114 = 0.1e1 / t122;
	t119 = 0.1e1 / t125;
	t144 = t104 * t114 * t119;
	t153 = t122 * t125;
	t94 = atan2(-t104, -t153);
	t92 = sin(t94);
	t93 = cos(t94);
	t101 = t104 ^ 2;
	t115 = 0.1e1 / t122 ^ 2;
	t120 = 0.1e1 / t125 ^ 2;
	t99 = t101 * t115 * t120 + 0.1e1;
	t95 = 0.1e1 / t99;
	t166 = (t93 * t144 - t92) * t95 + t92;
	t87 = -t92 * t104 - t93 * t153;
	t84 = 0.1e1 / t87;
	t116 = 0.1e1 / t124;
	t117 = 0.1e1 / t124 ^ 2;
	t85 = 0.1e1 / t87 ^ 2;
	t149 = qJD(2) * t123;
	t158 = t125 * t92;
	t163 = t104 * t93;
	t143 = t120 * t149;
	t160 = t114 * t95;
	t134 = -t123 * t141 - t124 * t125;
	t142 = t124 * t156;
	t135 = -t126 * t123 - t125 * t142;
	t90 = -t135 * qJD(1) - t134 * qJD(2);
	t77 = (t104 * t143 + t119 * t90) * t160;
	t74 = -t77 * t163 - t92 * t90 + (t93 * t149 + t77 * t158) * t122;
	t165 = t74 * t84 * t85;
	t154 = t120 * t123;
	t136 = t104 * t154 - t119 * t134;
	t78 = t136 * t160;
	t164 = t77 * t78;
	t162 = t135 * t85;
	t161 = t135 * t93;
	t159 = t119 * t95;
	t157 = t92 * t135;
	t155 = t117 * t126;
	t151 = t126 * t125;
	t150 = qJD(1) * t126;
	t102 = t135 ^ 2;
	t81 = t102 * t85 + 0.1e1;
	t138 = qJD(2) * t156 + qJD(1);
	t88 = -qJD(1) * t139 - qJD(2) * t151 + t138 * t152;
	t148 = 0.2e1 * (-t102 * t165 + t88 * t162) / t81 ^ 2;
	t147 = 0.2e1 * t165;
	t121 = t119 * t120;
	t146 = -0.2e1 * (t101 * t121 * t149 + t104 * t120 * t90) * t115 / t99 ^ 2;
	t140 = t123 * t142;
	t108 = -t140 + t151;
	t103 = t108 ^ 2;
	t100 = t103 * t117 * t115 + 0.1e1;
	t118 = t116 * t117;
	t89 = t134 * qJD(1) + t135 * qJD(2);
	t145 = 0.2e1 * (-t103 * t118 * t150 + t108 * t117 * t89) * t115 / t100 ^ 2;
	t133 = t119 * t146 + t95 * t143;
	t97 = 0.1e1 / t100;
	t91 = -qJD(1) * t140 - t124 * t149 + t138 * t151;
	t79 = 0.1e1 / t81;
	t76 = t166 * t135;
	t75 = -t78 * t163 + t92 * t134 + (t123 * t93 + t78 * t158) * t122;
	t73 = (t136 * t146 + (t90 * t154 + t119 * t91 + (-t134 * t154 + (0.2e1 * t121 * t123 ^ 2 + t119) * t104) * qJD(2)) * t95) * t114;
	t1 = [(-t133 * t135 - t88 * t159) * t114, t73, 0, 0, 0, 0; t104 * t84 * t148 + (-t90 * t84 + (t104 * t74 + t76 * t88) * t85) * t79 - (t76 * t147 * t79 + (t76 * t148 + ((t77 * t95 * t144 + t146) * t157 + ((t95 - 0.1e1) * t77 + (-t133 * t104 - t90 * t159) * t114) * t161 - t166 * t88) * t79) * t85) * t135, (-t108 * t84 - t75 * t162) * t148 + (-t75 * t135 * t147 + t89 * t84 + (-t108 * t74 + t75 * t88 + (t104 * t164 - t91) * t157 + (-t104 * t73 + t134 * t77 - t78 * t90) * t161) * t85 + ((-qJD(2) * t78 - t77) * t92 * t123 + (t73 * t92 + (qJD(2) + t164) * t93) * t125) * t122 * t162) * t79, 0, 0, 0, 0; ((t108 * t155 - t116 * t134) * t145 + (-t89 * t155 - t116 * t91 + (-t134 * t155 + (0.2e1 * t118 * t126 ^ 2 + t116) * t108) * qJD(1)) * t97) * t114, (t116 * t88 * t97 - (t117 * t97 * t150 + t116 * t145) * t135) * t114, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:32
	% DurationCPUTime: 0.88s
	% Computational Cost: add. (1334->90), mult. (4303->199), div. (668->14), fcn. (5516->11), ass. (0->91)
	t174 = sin(qJ(2));
	t175 = sin(qJ(1));
	t177 = cos(qJ(2));
	t178 = cos(qJ(1));
	t225 = cos(pkin(6));
	t194 = t178 * t225;
	t157 = t174 * t194 + t175 * t177;
	t172 = sin(pkin(6));
	t213 = t172 * t174;
	t151 = atan2(-t157, t213);
	t147 = sin(t151);
	t148 = cos(t151);
	t154 = t157 ^ 2;
	t168 = 0.1e1 / t172 ^ 2;
	t170 = 0.1e1 / t174 ^ 2;
	t152 = t154 * t168 * t170 + 0.1e1;
	t149 = 0.1e1 / t152;
	t167 = 0.1e1 / t172;
	t169 = 0.1e1 / t174;
	t199 = t157 * t167 * t169;
	t226 = (t148 * t199 + t147) * t149 - t147;
	t131 = -t147 * t157 + t148 * t213;
	t128 = 0.1e1 / t131;
	t195 = t175 * t225;
	t159 = t178 * t174 + t177 * t195;
	t173 = sin(qJ(4));
	t176 = cos(qJ(4));
	t212 = t172 * t175;
	t144 = t159 * t173 + t176 * t212;
	t140 = 0.1e1 / t144;
	t129 = 0.1e1 / t131 ^ 2;
	t141 = 0.1e1 / t144 ^ 2;
	t190 = qJD(2) * t225 + qJD(1);
	t192 = t174 * t195;
	t207 = qJD(2) * t174;
	t209 = t178 * t177;
	t138 = -qJD(1) * t192 - t175 * t207 + t190 * t209;
	t206 = qJD(2) * t177;
	t196 = t170 * t206;
	t185 = -t138 * t169 + t157 * t196;
	t215 = t149 * t167;
	t120 = t185 * t215;
	t187 = -t147 * t213 - t148 * t157;
	t200 = t148 * t172 * t177;
	t116 = qJD(2) * t200 + t187 * t120 - t147 * t138;
	t224 = t116 * t128 * t129;
	t191 = t177 * t194;
	t210 = t175 * t174;
	t156 = -t191 + t210;
	t214 = t170 * t177;
	t186 = t156 * t169 + t157 * t214;
	t121 = t186 * t215;
	t117 = t187 * t121 + t147 * t156 + t200;
	t160 = -t192 + t209;
	t223 = t117 * t160;
	t135 = -qJD(1) * t191 - t178 * t206 + t190 * t210;
	t208 = qJD(1) * t172;
	t197 = t178 * t208;
	t126 = t144 * qJD(4) + t135 * t176 + t173 * t197;
	t143 = -t159 * t176 + t173 * t212;
	t139 = t143 ^ 2;
	t134 = t139 * t141 + 0.1e1;
	t218 = t141 * t143;
	t205 = qJD(4) * t143;
	t127 = -t135 * t173 + t176 * t197 - t205;
	t220 = t127 * t140 * t141;
	t222 = (t126 * t218 - t139 * t220) / t134 ^ 2;
	t171 = t169 * t170;
	t221 = (t138 * t157 * t170 - t154 * t171 * t206) * t168 / t152 ^ 2;
	t136 = t157 * qJD(1) + t159 * qJD(2);
	t219 = t136 * t129;
	t217 = t147 * t160;
	t216 = t148 * t160;
	t211 = t172 * t178;
	t155 = t160 ^ 2;
	t124 = t155 * t129 + 0.1e1;
	t204 = 0.2e1 * (-t155 * t224 - t160 * t219) / t124 ^ 2;
	t203 = 0.2e1 * t224;
	t202 = 0.2e1 * t222;
	t201 = -0.2e1 * t221;
	t198 = t175 * t208;
	t193 = 0.2e1 * t143 * t220;
	t188 = t176 * t140 + t173 * t218;
	t146 = -t156 * t173 + t176 * t211;
	t145 = t156 * t176 + t173 * t211;
	t137 = t159 * qJD(1) + qJD(2) * t157;
	t132 = 0.1e1 / t134;
	t122 = 0.1e1 / t124;
	t119 = t226 * t160;
	t115 = (t186 * t201 + (t138 * t214 + t137 * t169 + (-t156 * t214 + (-0.2e1 * t171 * t177 ^ 2 - t169) * t157) * qJD(2)) * t149) * t167;
	t1 = [(0.2e1 * t160 * t169 * t221 + (t136 * t169 + t160 * t196) * t149) * t167, t115, 0, 0, 0, 0; t157 * t128 * t204 + (-t138 * t128 + (t116 * t157 + t119 * t136) * t129) * t122 + ((t119 * t203 + t226 * t219) * t122 + (t119 * t204 + (-(-t120 * t149 * t199 + t201) * t217 - (t199 * t201 - t120 + (-t167 * t185 + t120) * t149) * t216) * t122) * t129) * t160, (t128 * t159 + t129 * t223) * t204 + (t203 * t223 + t135 * t128 + (t159 * t116 + t117 * t136 - (-t172 * t207 - t115 * t157 - t121 * t138 + (-t121 * t213 + t156) * t120) * t216 - (t120 * t121 * t157 + t137 + (-t115 * t174 + (-qJD(2) * t121 - t120) * t177) * t172) * t217) * t129) * t122, 0, 0, 0, 0; (-t140 * t145 + t146 * t218) * t202 + ((t146 * qJD(4) + t137 * t176 - t173 * t198) * t140 + t146 * t193 + (-t145 * t127 - (-t145 * qJD(4) - t137 * t173 - t176 * t198) * t143 - t146 * t126) * t141) * t132, t188 * t160 * t202 + (t188 * t136 + ((qJD(4) * t140 + t193) * t173 + (-t126 * t173 + (t127 - t205) * t176) * t141) * t160) * t132, 0, -0.2e1 * t222 + 0.2e1 * (t126 * t141 * t132 + (-t132 * t220 - t141 * t222) * t143) * t143, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:32
	% DurationCPUTime: 0.88s
	% Computational Cost: add. (1643->91), mult. (4303->199), div. (668->14), fcn. (5516->11), ass. (0->92)
	t184 = sin(qJ(2));
	t185 = sin(qJ(1));
	t186 = cos(qJ(2));
	t187 = cos(qJ(1));
	t234 = cos(pkin(6));
	t203 = t187 * t234;
	t165 = t184 * t203 + t185 * t186;
	t183 = sin(pkin(6));
	t222 = t183 * t184;
	t159 = atan2(-t165, t222);
	t155 = sin(t159);
	t156 = cos(t159);
	t162 = t165 ^ 2;
	t178 = 0.1e1 / t183 ^ 2;
	t181 = 0.1e1 / t184 ^ 2;
	t160 = t162 * t178 * t181 + 0.1e1;
	t157 = 0.1e1 / t160;
	t177 = 0.1e1 / t183;
	t180 = 0.1e1 / t184;
	t208 = t165 * t177 * t180;
	t235 = (t156 * t208 + t155) * t157 - t155;
	t139 = -t155 * t165 + t156 * t222;
	t136 = 0.1e1 / t139;
	t204 = t185 * t234;
	t167 = t187 * t184 + t186 * t204;
	t179 = qJ(4) + pkin(11);
	t175 = sin(t179);
	t176 = cos(t179);
	t221 = t183 * t185;
	t152 = t167 * t175 + t176 * t221;
	t148 = 0.1e1 / t152;
	t137 = 0.1e1 / t139 ^ 2;
	t149 = 0.1e1 / t152 ^ 2;
	t199 = qJD(2) * t234 + qJD(1);
	t201 = t184 * t204;
	t216 = qJD(2) * t184;
	t218 = t187 * t186;
	t146 = -qJD(1) * t201 - t185 * t216 + t199 * t218;
	t215 = qJD(2) * t186;
	t205 = t181 * t215;
	t194 = -t146 * t180 + t165 * t205;
	t224 = t157 * t177;
	t128 = t194 * t224;
	t196 = -t155 * t222 - t156 * t165;
	t209 = t156 * t183 * t186;
	t124 = qJD(2) * t209 + t196 * t128 - t155 * t146;
	t233 = t124 * t136 * t137;
	t200 = t186 * t203;
	t219 = t185 * t184;
	t164 = -t200 + t219;
	t223 = t181 * t186;
	t195 = t164 * t180 + t165 * t223;
	t129 = t195 * t224;
	t125 = t196 * t129 + t155 * t164 + t209;
	t168 = -t201 + t218;
	t232 = t125 * t168;
	t143 = -qJD(1) * t200 - t187 * t215 + t199 * t219;
	t217 = qJD(1) * t183;
	t206 = t187 * t217;
	t133 = t152 * qJD(4) + t143 * t176 + t175 * t206;
	t151 = -t167 * t176 + t175 * t221;
	t147 = t151 ^ 2;
	t142 = t147 * t149 + 0.1e1;
	t227 = t149 * t151;
	t214 = qJD(4) * t151;
	t134 = -t143 * t175 + t176 * t206 - t214;
	t230 = t134 * t148 * t149;
	t231 = (t133 * t227 - t147 * t230) / t142 ^ 2;
	t182 = t180 * t181;
	t229 = (t146 * t165 * t181 - t162 * t182 * t215) * t178 / t160 ^ 2;
	t144 = qJD(1) * t165 + qJD(2) * t167;
	t228 = t144 * t137;
	t226 = t155 * t168;
	t225 = t156 * t168;
	t220 = t183 * t187;
	t163 = t168 ^ 2;
	t132 = t163 * t137 + 0.1e1;
	t213 = 0.2e1 * (-t163 * t233 - t168 * t228) / t132 ^ 2;
	t212 = 0.2e1 * t233;
	t211 = 0.2e1 * t231;
	t210 = -0.2e1 * t229;
	t207 = t185 * t217;
	t202 = 0.2e1 * t151 * t230;
	t197 = t176 * t148 + t175 * t227;
	t154 = -t164 * t175 + t176 * t220;
	t153 = t164 * t176 + t175 * t220;
	t145 = qJD(1) * t167 + qJD(2) * t165;
	t140 = 0.1e1 / t142;
	t130 = 0.1e1 / t132;
	t127 = t235 * t168;
	t123 = (t195 * t210 + (t146 * t223 + t145 * t180 + (-t164 * t223 + (-0.2e1 * t182 * t186 ^ 2 - t180) * t165) * qJD(2)) * t157) * t177;
	t1 = [(0.2e1 * t168 * t180 * t229 + (t144 * t180 + t168 * t205) * t157) * t177, t123, 0, 0, 0, 0; t165 * t136 * t213 + (-t146 * t136 + (t124 * t165 + t127 * t144) * t137) * t130 + ((t127 * t212 + t235 * t228) * t130 + (t127 * t213 + (-(-t128 * t157 * t208 + t210) * t226 - (t208 * t210 - t128 + (-t177 * t194 + t128) * t157) * t225) * t130) * t137) * t168, (t136 * t167 + t137 * t232) * t213 + (t212 * t232 + t143 * t136 + (t167 * t124 + t125 * t144 - (-t183 * t216 - t123 * t165 - t129 * t146 + (-t129 * t222 + t164) * t128) * t225 - (t128 * t129 * t165 + t145 + (-t123 * t184 + (-qJD(2) * t129 - t128) * t186) * t183) * t226) * t137) * t130, 0, 0, 0, 0; (-t148 * t153 + t154 * t227) * t211 + ((t154 * qJD(4) + t145 * t176 - t175 * t207) * t148 + t154 * t202 + (-t153 * t134 - (-t153 * qJD(4) - t145 * t175 - t176 * t207) * t151 - t154 * t133) * t149) * t140, t197 * t168 * t211 + (t197 * t144 + ((qJD(4) * t148 + t202) * t175 + (-t133 * t175 + (t134 - t214) * t176) * t149) * t168) * t140, 0, -0.2e1 * t231 + 0.2e1 * (t133 * t149 * t140 + (-t140 * t230 - t149 * t231) * t151) * t151, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:34
	% DurationCPUTime: 2.59s
	% Computational Cost: add. (8428->153), mult. (13478->302), div. (726->12), fcn. (17045->13), ass. (0->129)
	t267 = sin(qJ(2));
	t268 = sin(qJ(1));
	t270 = cos(qJ(2));
	t271 = cos(qJ(1));
	t349 = cos(pkin(6));
	t301 = t271 * t349;
	t253 = t267 * t268 - t270 * t301;
	t264 = qJ(4) + pkin(11);
	t262 = sin(t264);
	t263 = cos(t264);
	t265 = sin(pkin(6));
	t327 = t265 * t271;
	t290 = t253 * t263 + t262 * t327;
	t238 = t290 ^ 2;
	t328 = t265 * t270;
	t287 = -t349 * t262 - t263 * t328;
	t248 = 0.1e1 / t287 ^ 2;
	t227 = t238 * t248 + 0.1e1;
	t221 = 0.1e1 / t227;
	t324 = qJD(1) * t268;
	t308 = t265 * t324;
	t286 = -t267 * t301 - t268 * t270;
	t302 = t268 * t349;
	t288 = t271 * t267 + t270 * t302;
	t279 = t288 * qJD(1) - t286 * qJD(2);
	t309 = t263 * t327;
	t357 = -qJD(4) * t309 - t279 * t263;
	t210 = (qJD(4) * t253 + t308) * t262 + t357;
	t251 = -t262 * t328 + t349 * t263;
	t323 = qJD(2) * t267;
	t306 = t265 * t323;
	t236 = t251 * qJD(4) - t263 * t306;
	t247 = 0.1e1 / t287;
	t335 = t290 * t248;
	t294 = t210 * t247 - t236 * t335;
	t194 = t294 * t221;
	t228 = atan2(t290, -t287);
	t215 = sin(t228);
	t216 = cos(t228);
	t296 = t215 * t287 + t216 * t290;
	t189 = t296 * t194 - t210 * t215 + t216 * t236;
	t206 = t215 * t290 - t216 * t287;
	t204 = 0.1e1 / t206 ^ 2;
	t358 = t189 * t204;
	t321 = qJD(4) * t262;
	t303 = t265 * t321;
	t320 = qJD(4) * t263;
	t211 = t253 * t320 + t279 * t262 + t263 * t308 + t271 * t303;
	t329 = t265 * t268;
	t240 = t288 * t262 + t263 * t329;
	t297 = t267 * t302;
	t326 = t271 * t270;
	t255 = -t297 + t326;
	t266 = sin(qJ(6));
	t269 = cos(qJ(6));
	t223 = t240 * t266 - t255 * t269;
	t356 = 0.2e1 * t223;
	t203 = 0.1e1 / t206;
	t355 = t203 * t358;
	t284 = t288 * t263;
	t239 = t262 * t329 - t284;
	t298 = 0.2e1 * t239 * t355;
	t282 = (t349 * qJD(1) + qJD(2)) * t326 - qJD(2) * t297 - t267 * t324;
	t307 = qJD(1) * t327;
	t213 = t240 * qJD(4) + t262 * t307 - t282 * t263;
	t341 = t213 * t204;
	t354 = -t341 + t298;
	t353 = t236 * t248;
	t330 = t265 * t267;
	t310 = t290 * t330;
	t289 = t247 * t286 + t248 * t310;
	t352 = t263 * t289;
	t351 = -t255 * t262 * qJD(6) - t282;
	t350 = -qJD(6) * t288 + t255 * t320;
	t224 = t240 * t269 + t255 * t266;
	t218 = 0.1e1 / t224;
	t219 = 0.1e1 / t224 ^ 2;
	t237 = t239 ^ 2;
	t200 = t204 * t237 + 0.1e1;
	t348 = (-t237 * t355 + t239 * t341) / t200 ^ 2;
	t214 = qJD(4) * t284 + t282 * t262 + t263 * t307 - t268 * t303;
	t233 = t286 * qJD(1) - t288 * qJD(2);
	t201 = t224 * qJD(6) + t214 * t266 - t233 * t269;
	t217 = t223 ^ 2;
	t209 = t217 * t219 + 0.1e1;
	t338 = t219 * t223;
	t319 = qJD(6) * t223;
	t202 = t214 * t269 + t233 * t266 - t319;
	t344 = t202 * t218 * t219;
	t347 = (t201 * t338 - t217 * t344) / t209 ^ 2;
	t337 = t247 * t353;
	t345 = (-t210 * t335 + t238 * t337) / t227 ^ 2;
	t343 = t204 * t239;
	t207 = 0.1e1 / t209;
	t342 = t207 * t219;
	t340 = t215 * t239;
	t339 = t216 * t239;
	t336 = t290 * t247;
	t334 = t290 * t251;
	t333 = t255 * t263;
	t332 = t262 * t266;
	t331 = t262 * t269;
	t322 = qJD(2) * t270;
	t318 = 0.2e1 * t348;
	t317 = -0.2e1 * t347;
	t316 = -0.2e1 * t345;
	t315 = 0.2e1 * t345;
	t313 = t219 * t347;
	t312 = t201 * t342;
	t311 = t223 * t344;
	t300 = t247 * t315;
	t299 = 0.2e1 * t311;
	t291 = -t253 * t262 + t309;
	t226 = t266 * t286 + t269 * t291;
	t225 = t266 * t291 - t269 * t286;
	t293 = -t218 * t266 + t269 * t338;
	t292 = t247 * t291 + t248 * t334;
	t285 = -t215 + (t216 * t336 + t215) * t221;
	t235 = t287 * qJD(4) + t262 * t306;
	t234 = -qJD(1) * t297 - t268 * t323 + (qJD(2) * t349 + qJD(1)) * t326;
	t230 = t255 * t331 - t288 * t266;
	t198 = 0.1e1 / t200;
	t197 = t221 * t352;
	t195 = t292 * t221;
	t191 = (-t215 * t286 - t216 * t330) * t263 + t296 * t197;
	t190 = -t296 * t195 + t215 * t291 + t216 * t251;
	t188 = t292 * t315 + (-0.2e1 * t334 * t337 + t211 * t247 + (t210 * t251 - t235 * t290 - t236 * t291) * t248) * t221;
	t186 = t316 * t352 + (-t289 * t321 + (0.2e1 * t310 * t337 - t234 * t247 + (t236 * t286 + (-t210 * t267 + t290 * t322) * t265) * t248) * t263) * t221;
	t1 = [-t239 * t300 + (t213 * t247 + t239 * t353) * t221, t186, 0, t188, 0, 0; -0.2e1 * t290 * t203 * t348 + ((-t253 * t321 - t262 * t308 - t357) * t203 - t290 * t358 - (t285 * t213 + ((-t194 * t221 * t336 + t316) * t215 + (-t290 * t300 - t194 + (t194 - t294) * t221) * t216) * t239) * t343) * t198 + (t354 * t198 + t343 * t318) * t285 * t239, (t191 * t343 + t203 * t333) * t318 + ((-t233 * t263 + t255 * t321) * t203 + t354 * t191 + (t333 * t189 - (t186 * t290 - t197 * t210 + (-t263 * t322 + t267 * t321) * t265 + (t197 * t287 - t263 * t286) * t194) * t339 - (t286 * t321 + t186 * t287 - t197 * t236 + t234 * t263 + (-t197 * t290 + t263 * t330) * t194) * t340) * t204) * t198, 0, (t190 * t343 - t203 * t240) * t318 + (t190 * t298 + t214 * t203 + (-t240 * t189 - t190 * t213 - (t188 * t290 + t195 * t210 + t235 + (-t195 * t287 + t291) * t194) * t339 - (t188 * t287 + t195 * t236 - t211 + (t195 * t290 - t251) * t194) * t340) * t204) * t198, 0, 0; 0.2e1 * (-t218 * t225 + t226 * t338) * t347 + ((t226 * qJD(6) - t211 * t266 + t234 * t269) * t218 + t226 * t299 + (-t225 * t202 - (-t225 * qJD(6) - t211 * t269 - t234 * t266) * t223 - t226 * t201) * t219) * t207, (t313 * t356 - t312) * t230 + (-t202 * t342 + t218 * t317) * (t255 * t332 + t288 * t269) + (t230 * t299 + (t332 * t218 - t331 * t338) * t233 + (-t351 * t218 - t350 * t338) * t269 + (t350 * t218 - t351 * t338) * t266) * t207, 0, t293 * t239 * t317 + (t293 * t213 + ((-qJD(6) * t218 - 0.2e1 * t311) * t269 + (t201 * t269 + (t202 - t319) * t266) * t219) * t239) * t207, 0, t317 + (t312 + (-t207 * t344 - t313) * t223) * t356;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end