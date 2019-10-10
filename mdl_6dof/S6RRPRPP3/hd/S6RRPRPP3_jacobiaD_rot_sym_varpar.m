% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
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
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RRPRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:59
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:11
	% DurationCPUTime: 1.03s
	% Computational Cost: add. (892->82), mult. (2191->191), div. (456->12), fcn. (2616->9), ass. (0->85)
	t100 = sin(qJ(2));
	t92 = t100 ^ 2;
	t102 = cos(qJ(2));
	t95 = 0.1e1 / t102 ^ 2;
	t141 = t92 * t95;
	t101 = sin(qJ(1));
	t122 = 0.1e1 + t141;
	t93 = t101 ^ 2;
	t90 = t93 * t141 + 0.1e1;
	t88 = 0.1e1 / t90;
	t113 = t122 * t88;
	t74 = t101 * t113;
	t157 = t101 * t74 - 0.1e1;
	t103 = cos(qJ(1));
	t127 = qJD(2) * t103;
	t118 = t100 * t127;
	t134 = t101 * t102;
	t98 = sin(pkin(9));
	t99 = cos(pkin(9));
	t82 = t103 * t98 - t99 * t134;
	t76 = t82 * qJD(1) - t99 * t118;
	t133 = t102 * t103;
	t84 = t101 * t98 + t99 * t133;
	t78 = 0.1e1 / t84;
	t79 = 0.1e1 / t84 ^ 2;
	t80 = t78 * t79;
	t145 = t76 * t80;
	t81 = -t103 * t99 - t98 * t134;
	t75 = t81 * qJD(1) - t98 * t118;
	t146 = t75 * t79;
	t83 = -t101 * t99 + t98 * t133;
	t77 = t83 ^ 2;
	t72 = t77 * t79 + 0.1e1;
	t156 = (-t77 * t145 + t83 * t146) / t72 ^ 2;
	t155 = t100 * t141;
	t143 = t83 * t99;
	t112 = t79 * t143 - t78 * t98;
	t70 = 0.1e1 / t72;
	t154 = t112 * t70;
	t135 = t101 * t100;
	t87 = atan2(-t135, -t102);
	t85 = sin(t87);
	t123 = t85 * t135;
	t86 = cos(t87);
	t69 = -t102 * t86 - t123;
	t66 = 0.1e1 / t69;
	t94 = 0.1e1 / t102;
	t67 = 0.1e1 / t69 ^ 2;
	t153 = 0.2e1 * t100;
	t152 = t88 - 0.1e1;
	t130 = qJD(1) * t103;
	t114 = t101 * t92 * t130;
	t128 = qJD(2) * t102;
	t97 = t103 ^ 2;
	t140 = t92 * t97;
	t129 = qJD(2) * t101;
	t138 = t102 * t85;
	t62 = (-(-t100 * t130 - t101 * t128) * t94 + t129 * t141) * t88;
	t57 = (t62 - t129) * t138 + (-t85 * t130 + (-t101 * t62 + qJD(2)) * t86) * t100;
	t150 = t57 * t66 * t67;
	t65 = t67 * t140 + 0.1e1;
	t151 = (-t140 * t150 + (t100 * t97 * t128 - t114) * t67) / t65 ^ 2;
	t63 = 0.1e1 / t65;
	t148 = t63 * t67;
	t111 = qJD(2) * (t100 + t155) * t94;
	t147 = (t93 * t111 + t95 * t114) / t90 ^ 2;
	t144 = t82 * t83;
	t142 = t88 * t94;
	t137 = t103 * t67;
	t136 = qJD(2) * t74;
	t132 = qJD(1) * t100;
	t131 = qJD(1) * t101;
	t126 = 0.2e1 * t150;
	t125 = t66 * t151;
	t124 = t101 * t142;
	t121 = t100 * t152;
	t120 = t63 * t128;
	t119 = t100 * t129;
	t117 = 0.2e1 * t67 * t151;
	t116 = -0.2e1 * t94 * t147;
	t115 = t92 * t124;
	t61 = (-t86 * t115 + t85 * t121) * t103;
	t59 = (-t101 + t74) * t138 - t157 * t86 * t100;
	t58 = t113 * t130 + 0.2e1 * (t111 * t88 - t122 * t147) * t101;
	t1 = [-t124 * t132 + (qJD(2) * t113 + t100 * t116) * t103, t58, 0, 0, 0, 0; (-t66 * t120 + (0.2e1 * t125 + (qJD(1) * t61 + t57) * t148) * t100) * t101 + (-t61 * t67 * t120 + (t61 * t117 + (t61 * t126 + ((-t62 * t115 - t152 * t128 + t147 * t153) * t85 + (-t62 * t121 + (t92 * t116 + (t153 + t155) * t88 * qJD(2)) * t101) * t86) * t137) * t63) * t100 + (-t66 + (-(t93 - t97) * t92 * t86 * t142 + t152 * t123) * t67) * t63 * t132) * t103, (-t66 * t63 * t131 + (-0.2e1 * t125 + (-qJD(2) * t59 - t57) * t148) * t103) * t102 + (t59 * t103 * t117 + (-t66 * t127 - ((-t101 * t58 - t130 * t74) * t86 + (t157 * t62 + t129 - t136) * t85) * t100 * t137 + (t103 * t126 + t67 * t131) * t59) * t63 - ((t58 - t130) * t85 + (t62 * t74 + qJD(2) + (-t62 - t136) * t101) * t86) * t133 * t148) * t100, 0, 0, 0, 0; 0.2e1 * (t79 * t144 - t78 * t81) * t156 + ((-t83 * qJD(1) + t98 * t119) * t78 + 0.2e1 * t144 * t145 + (-t81 * t76 - (-t84 * qJD(1) + t99 * t119) * t83 - t82 * t75) * t79) * t70, t102 * t127 * t154 + (-t131 * t154 + (-0.2e1 * t112 * t156 + (t99 * t146 + (-0.2e1 * t80 * t143 + t79 * t98) * t76) * t70) * t103) * t100, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:11
	% DurationCPUTime: 1.02s
	% Computational Cost: add. (1350->93), mult. (2519->209), div. (480->12), fcn. (2968->9), ass. (0->94)
	t133 = sin(qJ(1));
	t127 = t133 ^ 2;
	t132 = sin(qJ(2));
	t126 = t132 ^ 2;
	t134 = cos(qJ(2));
	t129 = 0.1e1 / t134 ^ 2;
	t182 = t126 * t129;
	t119 = t127 * t182 + 0.1e1;
	t125 = t132 * t126;
	t128 = 0.1e1 / t134;
	t179 = t128 * t132;
	t144 = qJD(2) * (t125 * t128 * t129 + t179);
	t135 = cos(qJ(1));
	t171 = qJD(1) * t135;
	t180 = t126 * t133;
	t148 = t171 * t180;
	t188 = (t127 * t144 + t129 * t148) / t119 ^ 2;
	t198 = -0.2e1 * t188;
	t124 = pkin(9) + qJ(4);
	t123 = cos(t124);
	t173 = t134 * t135;
	t122 = sin(t124);
	t177 = t133 * t122;
	t112 = t123 * t173 + t177;
	t107 = 0.1e1 / t112 ^ 2;
	t176 = t133 * t123;
	t111 = t122 * t173 - t176;
	t185 = t111 * t123;
	t106 = 0.1e1 / t112;
	t187 = t106 * t122;
	t146 = t107 * t185 - t187;
	t105 = t111 ^ 2;
	t98 = t105 * t107 + 0.1e1;
	t96 = 0.1e1 / t98;
	t197 = t146 * t96;
	t155 = 0.1e1 + t182;
	t196 = t133 * t155;
	t175 = t133 * t132;
	t116 = atan2(-t175, -t134);
	t114 = cos(t116);
	t113 = sin(t116);
	t161 = t113 * t175;
	t102 = -t114 * t134 - t161;
	t99 = 0.1e1 / t102;
	t100 = 0.1e1 / t102 ^ 2;
	t117 = 0.1e1 / t119;
	t195 = t117 - 0.1e1;
	t131 = t135 ^ 2;
	t169 = qJD(2) * t134;
	t181 = t126 * t131;
	t170 = qJD(2) * t133;
	t157 = t129 * t170;
	t158 = t132 * t171;
	t92 = (-(-t133 * t169 - t158) * t128 + t126 * t157) * t117;
	t152 = t92 - t170;
	t153 = -t133 * t92 + qJD(2);
	t184 = t114 * t132;
	t86 = t153 * t184 + (t152 * t134 - t158) * t113;
	t191 = t99 * t100 * t86;
	t95 = t100 * t181 + 0.1e1;
	t194 = (-t181 * t191 + (t131 * t132 * t169 - t148) * t100) / t95 ^ 2;
	t186 = t107 * t111;
	t150 = -qJD(1) * t134 + qJD(4);
	t151 = qJD(4) * t134 - qJD(1);
	t168 = qJD(2) * t135;
	t156 = t132 * t168;
	t183 = t122 * t135;
	t91 = -t151 * t183 + (t150 * t133 - t156) * t123;
	t190 = t106 * t107 * t91;
	t174 = t133 * t134;
	t145 = t122 * t174 + t123 * t135;
	t90 = t145 * qJD(1) - qJD(4) * t112 + t122 * t156;
	t193 = (-t105 * t190 - t90 * t186) / t98 ^ 2;
	t93 = 0.1e1 / t95;
	t192 = t100 * t93;
	t178 = t132 * t135;
	t172 = qJD(1) * t133;
	t167 = 0.2e1 * t194;
	t166 = -0.2e1 * t193;
	t165 = 0.2e1 * t191;
	t164 = t99 * t194;
	t163 = t111 * t190;
	t162 = t93 * t169;
	t160 = t117 * t126 * t128;
	t154 = t128 * t198;
	t149 = t133 * t160;
	t147 = t155 * t135;
	t143 = t132 * t170 + t150 * t135;
	t110 = -t123 * t174 + t183;
	t104 = t117 * t196;
	t89 = (t195 * t132 * t113 - t114 * t149) * t135;
	t88 = -t113 * t174 + t184 + (t113 * t134 - t114 * t175) * t104;
	t87 = t196 * t198 + (qJD(1) * t147 + 0.2e1 * t133 * t144) * t117;
	t1 = [t154 * t178 + (qJD(2) * t147 - t172 * t179) * t117, t87, 0, 0, 0, 0; (-t99 * t162 + (0.2e1 * t164 + (qJD(1) * t89 + t86) * t192) * t132) * t133 + ((-t89 * t162 + (t89 * t167 + ((0.2e1 * t132 * t188 - t92 * t149 - t195 * t169) * t113 + (t154 * t180 + t132 * t92 + (t125 * t157 - (t92 - 0.2e1 * t170) * t132) * t117) * t114) * t93 * t135) * t132) * t100 + (t89 * t165 + (-t99 + ((-t127 + t131) * t114 * t160 + t195 * t161) * t100) * qJD(1)) * t132 * t93) * t135, (-t99 * t93 * t172 + (-0.2e1 * t164 + (-qJD(2) * t88 - t86) * t192) * t135) * t134 + (t88 * t135 * t100 * t167 + ((-qJD(2) * t99 + t88 * t165) * t135 + (t88 * t172 + (-(-t104 * t171 - t133 * t87) * t114 - ((t104 * t133 - 0.1e1) * t92 + (-t104 + t133) * qJD(2)) * t113) * t178) * t100) * t93 - ((t87 - t171) * t113 + (t152 * t104 + t153) * t114) * t173 * t192) * t132, 0, 0, 0, 0; 0.2e1 * (t106 * t145 + t110 * t186) * t193 + (0.2e1 * t110 * t163 - t151 * t106 * t176 + t143 * t187 + (-t151 * t111 * t177 + t110 * t90 - t143 * t185 + t145 * t91) * t107) * t96, t134 * t168 * t197 + (-t172 * t197 + (t146 * t166 + ((-qJD(4) * t106 - 0.2e1 * t163) * t123 + (-t123 * t90 + (-qJD(4) * t111 + t91) * t122) * t107) * t96) * t135) * t132, 0, t166 + 0.2e1 * (-t107 * t90 * t96 + (-t107 * t193 - t96 * t190) * t111) * t111, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:12
	% DurationCPUTime: 1.62s
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
	t1 = [t233 * t121 + 0.2e1 * t134 * t178, t98, 0, t96, 0, 0; t131 * t182 + (-t117 * t111 + (t102 * t115 + t131 * t99) * t112) * t107 + (t102 * t181 + (t102 * t201 + (t115 * t121 - t115 - (-t103 * t121 * t194 + t200) * t134) * t112 * t119 + (-(-0.2e1 * t131 * t178 - t103) * t230 + (-(t103 + t195) * t134 + t233 * t131) * t112 * t121) * t120) * t107) * t134, t101 * t134 * t181 + (-(-t98 * t225 + (t103 * t228 - t117 * t120) * t110) * t230 + t176 * t101 + (-t111 * t235 - (-t110 * t226 + t119 * t214 + t223) * t230) * t206) * t107 + (t182 * t235 + ((-t111 * t207 - (t211 * qJD(2) - t103) * t197) * t161 + (t111 * t210 + (t162 * t99 - (-t98 + t209) * t227 - (t211 * t103 - qJD(2)) * t224) * t112) * t159) * t107) * t149, 0, (t100 * t230 - t111 * t135) * t202 + (t116 * t111 + t176 * t100 - (-t118 + (-t103 * t150 - t149 * t96) * t159 - t171 * t106) * t197 + (-t135 * t99 - (t150 * t208 - t149 * t205 - t106 * t117 - t131 * t96 + (-t106 * t215 - t133) * t103) * t224) * t112) * t107, 0, 0; t174 * t180 + (t174 * t189 + (t116 * t216 - t118 * t156 + (-t133 * t216 + (0.2e1 * t158 * t160 ^ 2 + t156) * t135) * qJD(1)) * t152) * t125, (t135 * t191 + t150) * t199 + (-t116 * t191 + qJD(4) * t149 + (t156 * t234 - t190 * t217) * t135) * t125, 0, t134 * t156 * t180 + (t115 * t152 * t156 + (-t152 * t190 + t156 * t189) * t134) * t125, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:59:10
	% EndTime: 2019-10-10 09:59:12
	% DurationCPUTime: 1.42s
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
	t230 = t134 * (-t142 * t147 * t201 + t141 * t184) + t115 * t219;
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
	t1 = [t230 * t120 + 0.2e1 * t134 * t172, t97, 0, t95, 0, 0; t131 * t179 + (-t117 * t110 + (t101 * t115 + t131 * t98) * t111) * t106 + (t101 * t178 + (t101 * t196 + (t115 * t120 - t115 - (-t102 * t120 * t189 + t195) * t134) * t111 * t118 + (-(-0.2e1 * t131 * t172 - t102) * t227 + (-(t102 + t190) * t134 + t230 * t131) * t111 * t120) * t119) * t106) * t134, t100 * t134 * t178 + (-(-t97 * t222 + (t102 * t225 - t117 * t119) * t109) * t227 + (t177 + t226) * t100 + (t110 * t232 - (t109 * t223 - t118 * t211 - t220) * t227) * t201) * t106 + (t179 * t232 + ((-t110 * t202 - (qJD(2) * t207 - t102) * t192) * t156 + (t110 * t206 + (t157 * t98 - (-t97 + t204) * t224 - (t102 * t207 - qJD(2)) * t221) * t111) * t154) * t106) * t145, 0, (-t110 * t133 + t99 * t227) * t197 + (t99 * t177 + t114 * t110 - (t116 + (t102 * t144 - t145 * t95) * t154 + t166 * t105) * t192 + (t99 * t115 - t133 * t98 - (-t144 * t203 - t145 * t199 + t105 * t117 - t131 * t95 + (t105 * t212 + t130) * t102) * t221) * t111) * t106, 0, 0; t169 * t176 + (t169 * t184 + (t114 * t213 + t116 * t151 + (t130 * t213 + (0.2e1 * t153 * t155 ^ 2 + t151) * t133) * qJD(1)) * t147) * t124, (t133 * t186 - t144) * t194 + (-t114 * t186 + t200 + (t151 * t231 - t185 * t215) * t133) * t124, 0, t134 * t151 * t176 + (t115 * t147 * t151 + (-t147 * t152 * t206 + t151 * t184) * t134) * t124, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end