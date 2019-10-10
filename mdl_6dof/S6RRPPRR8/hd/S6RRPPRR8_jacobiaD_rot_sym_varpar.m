% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR8
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
%   Wie in S6RRPPRR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:48
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR8_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR8_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:32
	% DurationCPUTime: 0.96s
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
	t98 = sin(pkin(10));
	t99 = cos(pkin(10));
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
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:32
	% DurationCPUTime: 0.78s
	% Computational Cost: add. (926->76), mult. (3228->188), div. (613->15), fcn. (4191->9), ass. (0->85)
	t120 = cos(pkin(10));
	t172 = 0.2e1 * t120;
	t121 = sin(qJ(2));
	t115 = 0.1e1 / t121;
	t119 = sin(pkin(10));
	t124 = cos(qJ(1));
	t149 = t124 * t120;
	t122 = sin(qJ(1));
	t123 = cos(qJ(2));
	t152 = t122 * t123;
	t104 = t119 * t152 + t149;
	t112 = 0.1e1 / t119;
	t158 = t104 * t112;
	t140 = t115 * t158;
	t154 = t121 * t119;
	t94 = atan2(-t104, t154);
	t90 = sin(t94);
	t91 = cos(t94);
	t113 = 0.1e1 / t119 ^ 2;
	t116 = 0.1e1 / t121 ^ 2;
	t99 = t104 ^ 2;
	t97 = t113 * t116 * t99 + 0.1e1;
	t92 = 0.1e1 / t97;
	t171 = (t140 * t91 + t90) * t92 - t90;
	t86 = -t104 * t90 + t154 * t91;
	t83 = 0.1e1 / t86;
	t150 = t124 * t119;
	t138 = t123 * t150;
	t107 = -t120 * t122 + t138;
	t75 = t171 * t107;
	t170 = 0.2e1 * t75;
	t153 = t122 * t119;
	t108 = t123 * t149 + t153;
	t101 = 0.1e1 / t108;
	t102 = 0.1e1 / t108 ^ 2;
	t84 = 0.1e1 / t86 ^ 2;
	t100 = t107 ^ 2;
	t145 = qJD(2) * t124;
	t136 = t121 * t145;
	t87 = qJD(1) * t104 + t119 * t136;
	t166 = t84 * t87;
	t146 = qJD(2) * t123;
	t161 = t121 * t90;
	t164 = t104 * t91;
	t137 = t116 * t146;
	t148 = qJD(1) * t122;
	t89 = -qJD(2) * t121 * t153 + qJD(1) * t138 - t120 * t148;
	t76 = (t104 * t137 - t115 * t89) * t92 * t112;
	t73 = -t76 * t164 - t90 * t89 + (t146 * t91 - t161 * t76) * t119;
	t168 = t73 * t83 * t84;
	t80 = t100 * t84 + 0.1e1;
	t169 = (-t100 * t168 - t107 * t166) / t80 ^ 2;
	t114 = t121 ^ 2;
	t117 = t115 / t114;
	t167 = (t104 * t116 * t89 - t117 * t146 * t99) * t113 / t97 ^ 2;
	t106 = -t120 * t152 + t150;
	t133 = t120 * t136;
	t88 = qJD(1) * t106 - t133;
	t165 = t101 * t102 * t88;
	t163 = t107 * t91;
	t162 = t115 * t92;
	t160 = t90 * t107;
	t155 = t116 * t123;
	t82 = (t155 * t158 + t122) * t92;
	t159 = t122 - t82;
	t157 = t106 * t124;
	t118 = t124 ^ 2;
	t156 = t114 * t118;
	t151 = t123 * t124;
	t147 = qJD(1) * t124;
	t130 = t114 * t122 * t147 - t118 * t121 * t146;
	t134 = t156 * t165;
	t139 = t102 * t156;
	t98 = 0.1e1 + t139;
	t144 = 0.2e1 * (-t102 * t130 - t134) / t98 ^ 2;
	t143 = -0.2e1 * t167;
	t142 = t84 * t169;
	t141 = t84 * t160;
	t135 = 0.2e1 * t83 * t169;
	t131 = 0.2e1 * t115 * t167 + t137 * t92;
	t95 = 0.1e1 / t98;
	t78 = 0.1e1 / t80;
	t74 = -t82 * t164 + (t123 * t91 + t159 * t161) * t119;
	t72 = t122 * t143 + t92 * t147 + (t89 * t92 * t155 + (t143 * t155 + (-0.2e1 * t117 * t123 ^ 2 - t115) * t92 * qJD(2)) * t104) * t112;
	t1 = [(t107 * t131 + t162 * t87) * t112, t72, 0, 0, 0, 0; t104 * t135 + (-t89 * t83 + (t104 * t73 + t75 * t87) * t84) * t78 + (t142 * t170 + (t168 * t170 - (-t140 * t76 * t92 + t143) * t141 - ((t92 - 0.1e1) * t76 + (-t104 * t131 + t162 * t89) * t112) * t84 * t163 + t171 * t166) * t78) * t107, t74 * t78 * t166 + (-(-t82 * t91 * t89 + (t76 * t82 * t90 - t72 * t91) * t104) * t84 * t78 + 0.2e1 * (t168 * t78 + t142) * t74) * t107 + (t124 * t135 * t121 + ((-t83 * t145 - (qJD(2) * t159 - t76) * t141) * t123 + (t83 * t148 + (t124 * t73 - (-t72 + t147) * t160 - (t159 * t76 - qJD(2)) * t163) * t84) * t121) * t78) * t119, 0, 0, 0, 0; (-t101 * t122 - t102 * t157) * t121 * t144 + (-0.2e1 * t121 * t157 * t165 + (t121 * t147 + t122 * t146) * t101 + (((-t88 + t133) * t122 - t108 * t147) * t121 + (-t121 * t148 + t123 * t145) * t106) * t102) * t95, (t101 * t151 + t120 * t139) * t144 + (t134 * t172 + (t123 * t148 + t136) * t101 + (t130 * t172 + t151 * t88) * t102) * t95, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:33
	% DurationCPUTime: 1.18s
	% Computational Cost: add. (926->105), mult. (3345->233), div. (468->12), fcn. (4023->11), ass. (0->105)
	t155 = sin(qJ(2));
	t146 = t155 ^ 2;
	t158 = cos(qJ(2));
	t149 = 0.1e1 / t158 ^ 2;
	t205 = t146 * t149;
	t156 = sin(qJ(1));
	t147 = t156 ^ 2;
	t144 = t147 * t205 + 0.1e1;
	t148 = 0.1e1 / t158;
	t202 = t148 * t155;
	t224 = t155 * t205;
	t167 = qJD(2) * (t148 * t224 + t202);
	t159 = cos(qJ(1));
	t195 = qJD(1) * t159;
	t203 = t146 * t156;
	t173 = t195 * t203;
	t208 = (t147 * t167 + t149 * t173) / t144 ^ 2;
	t225 = -0.2e1 * t208;
	t177 = 0.1e1 + t205;
	t223 = t156 * t177;
	t154 = sin(qJ(5));
	t157 = cos(qJ(5));
	t191 = qJD(5) * t159;
	t196 = qJD(1) * t156;
	t222 = t154 * t196 - t157 * t191;
	t221 = t154 * t191 + t157 * t196;
	t153 = cos(pkin(10));
	t152 = sin(pkin(10));
	t198 = t159 * t152;
	t137 = -t156 * t153 + t158 * t198;
	t197 = t159 * t153;
	t138 = t156 * t152 + t158 * t197;
	t121 = t137 * t154 + t138 * t157;
	t115 = 0.1e1 / t121;
	t200 = t156 * t155;
	t143 = atan2(t200, t158);
	t140 = cos(t143);
	t139 = sin(t143);
	t186 = t139 * t200;
	t125 = t140 * t158 + t186;
	t122 = 0.1e1 / t125;
	t116 = 0.1e1 / t121 ^ 2;
	t123 = 0.1e1 / t125 ^ 2;
	t220 = 0.2e1 * t155;
	t141 = 0.1e1 / t144;
	t219 = t141 - 0.1e1;
	t199 = t156 * t158;
	t135 = -t152 * t199 - t197;
	t192 = qJD(2) * t159;
	t178 = t155 * t192;
	t128 = t135 * qJD(1) - t152 * t178;
	t136 = -t153 * t199 + t198;
	t129 = t136 * qJD(1) - t153 * t178;
	t104 = t121 * qJD(5) - t128 * t157 + t129 * t154;
	t170 = t137 * t157 - t138 * t154;
	t114 = t170 ^ 2;
	t108 = t114 * t116 + 0.1e1;
	t212 = t116 * t170;
	t105 = t170 * qJD(5) + t128 * t154 + t129 * t157;
	t117 = t115 * t116;
	t215 = t105 * t117;
	t218 = 0.1e1 / t108 ^ 2 * (-t104 * t212 - t114 * t215);
	t151 = t159 ^ 2;
	t204 = t146 * t151;
	t113 = t123 * t204 + 0.1e1;
	t193 = qJD(2) * t158;
	t182 = t155 * t195;
	t194 = qJD(2) * t156;
	t110 = ((t156 * t193 + t182) * t148 + t194 * t205) * t141;
	t206 = t140 * t155;
	t101 = (t110 * t156 - qJD(2)) * t206 + (t182 + (-t110 + t194) * t158) * t139;
	t216 = t101 * t122 * t123;
	t217 = (-t204 * t216 + (t151 * t155 * t193 - t173) * t123) / t113 ^ 2;
	t214 = t110 * t139;
	t213 = t110 * t155;
	t168 = -t152 * t154 - t153 * t157;
	t201 = t155 * t159;
	t133 = t168 * t201;
	t211 = t116 * t133;
	t210 = t123 * t155;
	t209 = t123 * t159;
	t127 = t141 * t223;
	t207 = t127 * t156;
	t190 = 0.2e1 * t218;
	t189 = -0.2e1 * t216;
	t188 = -0.2e1 * t117 * t170;
	t187 = t123 * t201;
	t185 = t141 * t146 * t148;
	t181 = t155 * t194;
	t176 = -0.2e1 * t155 * t217;
	t175 = t148 * t225;
	t174 = t156 * t185;
	t172 = t177 * t159;
	t171 = t135 * t157 - t136 * t154;
	t119 = t135 * t154 + t136 * t157;
	t169 = t152 * t157 - t153 * t154;
	t132 = t169 * t201;
	t131 = -t138 * qJD(1) + t153 * t181;
	t130 = -t137 * qJD(1) + t152 * t181;
	t111 = 0.1e1 / t113;
	t109 = (-t219 * t155 * t139 + t140 * t174) * t159;
	t106 = 0.1e1 / t108;
	t103 = t139 * t199 - t206 + (-t139 * t158 + t140 * t200) * t127;
	t102 = t223 * t225 + (qJD(1) * t172 + 0.2e1 * t156 * t167) * t141;
	t1 = [t175 * t201 + (qJD(2) * t172 - t196 * t202) * t141, t102, 0, 0, 0, 0; (t122 * t176 + (t122 * t193 + (-qJD(1) * t109 - t101) * t210) * t111) * t156 + (t123 * t176 * t109 + (((-t110 * t174 - t219 * t193 + t208 * t220) * t139 + (t175 * t203 + t213 + (-t213 + (t220 + t224) * t194) * t141) * t140) * t187 + (t123 * t193 + t155 * t189) * t109 + (t122 + ((-t147 + t151) * t140 * t185 + t219 * t186) * t123) * t155 * qJD(1)) * t111) * t159, 0.2e1 * (-t103 * t210 + t122 * t158) * t159 * t217 + ((t122 * t196 + (qJD(2) * t103 + t101) * t209) * t158 + (t122 * t192 + (t102 * t140 * t156 - t139 * t194 - t207 * t214 + t214 + (qJD(2) * t139 + t140 * t195) * t127) * t187 + (-t123 * t196 + t159 * t189) * t103 + ((-t102 + t195) * t139 + ((-0.1e1 + t207) * qJD(2) + (-t127 + t156) * t110) * t140) * t158 * t209) * t155) * t111, 0, 0, 0, 0; (t115 * t171 - t119 * t212) * t190 + ((t119 * qJD(5) - t130 * t157 + t131 * t154) * t115 + t119 * t105 * t188 + (t171 * t105 + (t171 * qJD(5) + t130 * t154 + t131 * t157) * t170 - t119 * t104) * t116) * t106, (-t115 * t132 - t170 * t211) * t190 + (-t104 * t211 + (-t132 * t116 + t133 * t188) * t105 + (t169 * t115 + t168 * t212) * t158 * t192 + ((t222 * t115 + t221 * t212) * t153 + (-t221 * t115 + t222 * t212) * t152) * t155) * t106, 0, 0, -0.2e1 * t218 - 0.2e1 * (t104 * t116 * t106 - (-t106 * t215 - t116 * t218) * t170) * t170, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:32
	% EndTime: 2019-10-10 09:48:33
	% DurationCPUTime: 1.19s
	% Computational Cost: add. (1540->107), mult. (3805->232), div. (486->12), fcn. (4569->11), ass. (0->107)
	t187 = sin(qJ(2));
	t178 = t187 ^ 2;
	t189 = cos(qJ(2));
	t181 = 0.1e1 / t189 ^ 2;
	t237 = t178 * t181;
	t188 = sin(qJ(1));
	t179 = t188 ^ 2;
	t173 = t179 * t237 + 0.1e1;
	t180 = 0.1e1 / t189;
	t234 = t180 * t187;
	t256 = t187 * t237;
	t198 = qJD(2) * (t180 * t256 + t234);
	t190 = cos(qJ(1));
	t227 = qJD(1) * t190;
	t235 = t178 * t188;
	t202 = t227 * t235;
	t240 = (t179 * t198 + t181 * t202) / t173 ^ 2;
	t257 = -0.2e1 * t240;
	t210 = 0.1e1 + t237;
	t255 = t188 * t210;
	t176 = qJD(5) + qJD(6);
	t186 = cos(pkin(10));
	t228 = qJD(1) * t188;
	t185 = sin(pkin(10));
	t230 = t190 * t185;
	t254 = -t176 * t230 + t186 * t228;
	t229 = t190 * t186;
	t253 = t176 * t229 + t185 * t228;
	t166 = -t188 * t186 + t189 * t230;
	t167 = t188 * t185 + t189 * t229;
	t184 = qJ(5) + qJ(6);
	t174 = sin(t184);
	t175 = cos(t184);
	t147 = t166 * t174 + t167 * t175;
	t141 = 0.1e1 / t147;
	t232 = t188 * t187;
	t172 = atan2(t232, t189);
	t169 = cos(t172);
	t168 = sin(t172);
	t219 = t168 * t232;
	t154 = t169 * t189 + t219;
	t151 = 0.1e1 / t154;
	t142 = 0.1e1 / t147 ^ 2;
	t152 = 0.1e1 / t154 ^ 2;
	t252 = 0.2e1 * t187;
	t170 = 0.1e1 / t173;
	t251 = t170 - 0.1e1;
	t231 = t188 * t189;
	t165 = -t186 * t231 + t230;
	t224 = qJD(2) * t190;
	t211 = t187 * t224;
	t206 = t165 * qJD(1) + t166 * t176 - t186 * t211;
	t164 = -t185 * t231 - t229;
	t207 = -t164 * qJD(1) + t167 * t176 + t185 * t211;
	t133 = t206 * t174 + t207 * t175;
	t146 = -t166 * t175 + t167 * t174;
	t140 = t146 ^ 2;
	t137 = t140 * t142 + 0.1e1;
	t244 = t142 * t146;
	t134 = -t207 * t174 + t206 * t175;
	t143 = t141 * t142;
	t247 = t134 * t143;
	t250 = (t133 * t244 - t140 * t247) / t137 ^ 2;
	t183 = t190 ^ 2;
	t236 = t178 * t183;
	t150 = t152 * t236 + 0.1e1;
	t225 = qJD(2) * t189;
	t213 = t187 * t227;
	t226 = qJD(2) * t188;
	t139 = ((t188 * t225 + t213) * t180 + t226 * t237) * t170;
	t238 = t169 * t187;
	t130 = (t139 * t188 - qJD(2)) * t238 + (t213 + (-t139 + t226) * t189) * t168;
	t248 = t130 * t151 * t152;
	t249 = (-t236 * t248 + (t183 * t187 * t225 - t202) * t152) / t150 ^ 2;
	t246 = t139 * t168;
	t245 = t139 * t187;
	t199 = -t174 * t185 - t175 * t186;
	t233 = t187 * t190;
	t162 = t199 * t233;
	t243 = t142 * t162;
	t242 = t152 * t187;
	t241 = t152 * t190;
	t156 = t170 * t255;
	t239 = t156 * t188;
	t223 = 0.2e1 * t250;
	t222 = -0.2e1 * t248;
	t221 = 0.2e1 * t143 * t146;
	t220 = t152 * t233;
	t218 = t170 * t178 * t180;
	t212 = t187 * t226;
	t209 = -0.2e1 * t187 * t249;
	t208 = t180 * t257;
	t205 = t166 * qJD(1) + t165 * t176 - t185 * t212;
	t204 = -t167 * qJD(1) + t164 * t176 + t186 * t212;
	t203 = t188 * t218;
	t201 = t210 * t190;
	t200 = -t174 * t186 + t175 * t185;
	t161 = t200 * t233;
	t148 = 0.1e1 / t150;
	t145 = t164 * t174 + t165 * t175;
	t144 = -t164 * t175 + t165 * t174;
	t138 = (-t251 * t187 * t168 + t169 * t203) * t190;
	t135 = 0.1e1 / t137;
	t132 = t168 * t231 - t238 + (-t168 * t189 + t169 * t232) * t156;
	t131 = t255 * t257 + (qJD(1) * t201 + 0.2e1 * t188 * t198) * t170;
	t127 = -0.2e1 * t250 + 0.2e1 * (t133 * t142 * t135 + (-t135 * t247 - t142 * t250) * t146) * t146;
	t1 = [t208 * t233 + (qJD(2) * t201 - t228 * t234) * t170, t131, 0, 0, 0, 0; (t151 * t209 + (t151 * t225 + (-qJD(1) * t138 - t130) * t242) * t148) * t188 + (t152 * t209 * t138 + (((-t139 * t203 - t251 * t225 + t240 * t252) * t168 + (t208 * t235 + t245 + (-t245 + (t252 + t256) * t226) * t170) * t169) * t220 + (t152 * t225 + t187 * t222) * t138 + (t151 + ((-t179 + t183) * t169 * t218 + t251 * t219) * t152) * t187 * qJD(1)) * t148) * t190, 0.2e1 * (-t132 * t242 + t151 * t189) * t190 * t249 + ((t151 * t228 + (qJD(2) * t132 + t130) * t241) * t189 + (t151 * t224 + (t131 * t169 * t188 - t168 * t226 - t239 * t246 + t246 + (qJD(2) * t168 + t169 * t227) * t156) * t220 + (-t152 * t228 + t190 * t222) * t132 + ((-t131 + t227) * t168 + ((-0.1e1 + t239) * qJD(2) + (-t156 + t188) * t139) * t169) * t189 * t241) * t187) * t148, 0, 0, 0, 0; (-t141 * t144 + t145 * t244) * t223 + ((t174 * t204 + t175 * t205) * t141 + t145 * t134 * t221 + (-t144 * t134 - (-t174 * t205 + t175 * t204) * t146 - t145 * t133) * t142) * t135, (-t141 * t161 + t146 * t243) * t223 + (-t133 * t243 + (-t161 * t142 + t162 * t221) * t134 + (t141 * t200 - t199 * t244) * t189 * t224 + ((-t253 * t141 - t254 * t244) * t175 + (t254 * t141 - t253 * t244) * t174) * t187) * t135, 0, 0, t127, t127;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end