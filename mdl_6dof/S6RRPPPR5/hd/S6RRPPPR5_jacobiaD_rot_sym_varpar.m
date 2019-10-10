% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR5
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
%   Wie in S6RRPPPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:25
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:10
	% DurationCPUTime: 1.00s
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
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:10
	% DurationCPUTime: 0.78s
	% Computational Cost: add. (906->74), mult. (3124->186), div. (647->14), fcn. (4112->9), ass. (0->89)
	t119 = cos(qJ(2));
	t118 = sin(qJ(1));
	t159 = cos(pkin(9));
	t136 = t118 * t159;
	t116 = sin(pkin(9));
	t120 = cos(qJ(1));
	t154 = t120 * t116;
	t102 = t119 * t154 - t136;
	t135 = t120 * t159;
	t155 = t118 * t116;
	t103 = t119 * t135 + t155;
	t177 = qJD(1) * t103;
	t117 = sin(qJ(2));
	t109 = 0.1e1 / t117;
	t110 = 0.1e1 / t117 ^ 2;
	t111 = t109 * t110;
	t176 = qJD(2) * (0.2e1 * t111 * t119 ^ 2 + t109);
	t107 = 0.1e1 / t116;
	t99 = t119 * t155 + t135;
	t164 = t107 * t99;
	t145 = t109 * t164;
	t156 = t117 * t116;
	t91 = atan2(-t99, t156);
	t87 = sin(t91);
	t88 = cos(t91);
	t108 = 0.1e1 / t116 ^ 2;
	t96 = t99 ^ 2;
	t94 = t108 * t110 * t96 + 0.1e1;
	t89 = 0.1e1 / t94;
	t175 = (t88 * t145 + t87) * t89 - t87;
	t168 = t87 * t99;
	t83 = t88 * t156 - t168;
	t80 = 0.1e1 / t83;
	t72 = t175 * t102;
	t174 = 0.2e1 * t72;
	t113 = 0.1e1 / t120;
	t114 = 0.1e1 / t120 ^ 2;
	t81 = 0.1e1 / t83 ^ 2;
	t151 = qJD(2) * t117;
	t140 = t116 * t151;
	t84 = t99 * qJD(1) + t120 * t140;
	t169 = t81 * t84;
	t150 = qJD(2) * t119;
	t161 = t117 * t87;
	t167 = t88 * t99;
	t142 = t110 * t150;
	t130 = t99 * t142;
	t86 = t102 * qJD(1) - t118 * t140;
	t73 = (-t109 * t86 + t130) * t89 * t107;
	t70 = -t73 * t167 - t86 * t87 + (t88 * t150 - t73 * t161) * t116;
	t172 = t70 * t80 * t81;
	t97 = t102 ^ 2;
	t77 = t81 * t97 + 0.1e1;
	t173 = (-t102 * t169 - t97 * t172) / t77 ^ 2;
	t115 = t113 * t114;
	t141 = t111 * t150;
	t153 = qJD(1) * t118;
	t98 = t103 ^ 2;
	t162 = t114 * t98;
	t101 = -t119 * t136 + t154;
	t129 = t159 * t151;
	t85 = t101 * qJD(1) - t120 * t129;
	t95 = t110 * t162 + 0.1e1;
	t171 = (-t141 * t162 + (t103 * t114 * t85 + t115 * t98 * t153) * t110) / t95 ^ 2;
	t170 = (t110 * t86 * t99 - t96 * t141) * t108 / t94 ^ 2;
	t166 = t102 * t87;
	t165 = t102 * t88;
	t92 = 0.1e1 / t95;
	t163 = t113 * t92;
	t158 = t110 * t119;
	t79 = (t158 * t164 + t118) * t89;
	t160 = t118 - t79;
	t157 = t114 * t118;
	t152 = qJD(1) * t120;
	t149 = 0.2e1 * t171;
	t148 = -0.2e1 * t170;
	t147 = t81 * t173;
	t146 = t81 * t166;
	t144 = t92 * t157;
	t139 = 0.2e1 * t80 * t173;
	t138 = t99 * t148;
	t137 = -0.2e1 * t109 * t171;
	t133 = t109 * t144;
	t132 = t158 * t163;
	t128 = t144 * t158;
	t75 = 0.1e1 / t77;
	t71 = -t79 * t167 + (t119 * t88 + t160 * t161) * t116;
	t69 = t118 * t148 + t89 * t152 + (t138 * t158 + (t86 * t158 - t99 * t176) * t89) * t107;
	t1 = [(t109 * t84 * t89 + (0.2e1 * t109 * t170 + t89 * t142) * t102) * t107, t69, 0, 0, 0, 0; t99 * t139 + (-t86 * t80 + (t70 * t99 + t72 * t84) * t81) * t75 + (t147 * t174 + (t172 * t174 - (-t73 * t89 * t145 + t148) * t146 - ((t89 - 0.1e1) * t73 + (-t89 * t130 + (t86 * t89 + t138) * t109) * t107) * t81 * t165 + t175 * t169) * t75) * t102, t71 * t75 * t169 + (-(-t69 * t167 + (t73 * t168 - t86 * t88) * t79) * t81 * t75 + 0.2e1 * (t75 * t172 + t147) * t71) * t102 + (t120 * t139 * t117 + ((-t120 * qJD(2) * t80 - (t160 * qJD(2) - t73) * t146) * t119 + (t80 * t153 + (t120 * t70 - (-t69 + t152) * t166 - (t160 * t73 - qJD(2)) * t165) * t81) * t117) * t75) * t116, 0, 0, 0, 0; t85 * t133 + (qJD(1) * t133 - qJD(2) * t132 + t113 * t137) * t101 + ((t118 * t129 - t177) * t163 + (0.2e1 * t118 ^ 2 * t115 + t113) * t92 * t177) * t109 + (-qJD(2) * t128 + t137 * t157) * t103, -t85 * t132 + t159 * t149 + (-qJD(1) * t128 + (t149 * t158 + t92 * t176) * t113) * t103, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:10
	% DurationCPUTime: 0.78s
	% Computational Cost: add. (926->75), mult. (3228->188), div. (613->15), fcn. (4191->9), ass. (0->83)
	t122 = sin(qJ(2));
	t116 = 0.1e1 / t122;
	t121 = cos(pkin(9));
	t120 = sin(pkin(9));
	t125 = cos(qJ(1));
	t150 = t125 * t120;
	t123 = sin(qJ(1));
	t124 = cos(qJ(2));
	t152 = t123 * t124;
	t106 = t121 * t152 - t150;
	t113 = 0.1e1 / t121;
	t157 = t106 * t113;
	t140 = t116 * t157;
	t154 = t122 * t121;
	t95 = atan2(-t106, t154);
	t91 = sin(t95);
	t92 = cos(t95);
	t100 = t106 ^ 2;
	t114 = 0.1e1 / t121 ^ 2;
	t117 = 0.1e1 / t122 ^ 2;
	t98 = t100 * t114 * t117 + 0.1e1;
	t93 = 0.1e1 / t98;
	t171 = (t92 * t140 + t91) * t93 - t91;
	t87 = -t106 * t91 + t92 * t154;
	t84 = 0.1e1 / t87;
	t151 = t124 * t125;
	t109 = t123 * t120 + t121 * t151;
	t76 = t171 * t109;
	t170 = 0.2e1 * t76;
	t153 = t123 * t121;
	t108 = t124 * t150 - t153;
	t101 = 0.1e1 / t108;
	t102 = 0.1e1 / t108 ^ 2;
	t85 = 0.1e1 / t87 ^ 2;
	t104 = t109 ^ 2;
	t146 = qJD(2) * t125;
	t138 = t122 * t146;
	t149 = qJD(1) * t123;
	t133 = t124 * t149 + t138;
	t148 = qJD(1) * t125;
	t89 = -t120 * t148 + t133 * t121;
	t166 = t85 * t89;
	t147 = qJD(2) * t124;
	t160 = t122 * t91;
	t164 = t106 * t92;
	t139 = t117 * t147;
	t90 = -qJD(2) * t122 * t153 + t109 * qJD(1);
	t77 = (t106 * t139 - t116 * t90) * t93 * t113;
	t74 = -t77 * t164 - t90 * t91 + (t92 * t147 - t77 * t160) * t121;
	t168 = t74 * t84 * t85;
	t81 = t104 * t85 + 0.1e1;
	t169 = (-t104 * t168 - t109 * t166) / t81 ^ 2;
	t115 = t122 ^ 2;
	t118 = t116 / t115;
	t167 = (-t100 * t118 * t147 + t106 * t117 * t90) * t114 / t98 ^ 2;
	t105 = -t120 * t152 - t121 * t125;
	t135 = t120 * t138;
	t88 = t105 * qJD(1) - t135;
	t165 = t101 * t102 * t88;
	t163 = t109 * t91;
	t162 = t109 * t92;
	t161 = t116 * t93;
	t155 = t117 * t124;
	t83 = (t155 * t157 + t123) * t93;
	t159 = t123 - t83;
	t158 = t105 * t125;
	t119 = t125 ^ 2;
	t156 = t115 * t119;
	t131 = -t115 * t123 * t148 + t119 * t122 * t147;
	t136 = t156 * t165;
	t141 = t102 * t156;
	t99 = 0.1e1 + t141;
	t145 = 0.2e1 * (t131 * t102 - t136) / t99 ^ 2;
	t144 = -0.2e1 * t167;
	t143 = t85 * t169;
	t142 = t85 * t163;
	t137 = 0.2e1 * t84 * t169;
	t132 = 0.2e1 * t116 * t167 + t93 * t139;
	t96 = 0.1e1 / t99;
	t79 = 0.1e1 / t81;
	t75 = -t83 * t164 + (t124 * t92 + t159 * t160) * t121;
	t73 = t123 * t144 + t93 * t148 + (t90 * t93 * t155 + (t144 * t155 + (-0.2e1 * t118 * t124 ^ 2 - t116) * t93 * qJD(2)) * t106) * t113;
	t1 = [(t132 * t109 + t89 * t161) * t113, t73, 0, 0, 0, 0; t106 * t137 + (-t90 * t84 + (t106 * t74 + t76 * t89) * t85) * t79 + (t143 * t170 + (t168 * t170 - (-t77 * t93 * t140 + t144) * t142 - ((t93 - 0.1e1) * t77 + (-t132 * t106 + t90 * t161) * t113) * t85 * t162 + t171 * t166) * t79) * t109, t75 * t79 * t166 + (-(-t83 * t92 * t90 + (t77 * t83 * t91 - t73 * t92) * t106) * t85 * t79 + 0.2e1 * (t79 * t168 + t143) * t75) * t109 + (t125 * t137 * t122 + ((-t84 * t146 - (t159 * qJD(2) - t77) * t142) * t124 + (t84 * t149 + (t125 * t74 - (-t73 + t148) * t163 - (t159 * t77 - qJD(2)) * t162) * t85) * t122) * t79) * t121, 0, 0, 0, 0; (t101 * t123 + t102 * t158) * t122 * t145 + (0.2e1 * t122 * t158 * t165 + (-t122 * t148 - t123 * t147) * t101 + (((t88 - t135) * t123 + t108 * t148) * t122 + (t122 * t149 - t124 * t146) * t105) * t102) * t96, (-t101 * t151 - t120 * t141) * t145 + (-0.2e1 * t120 * t136 - t133 * t101 + (0.2e1 * t120 * t131 - t88 * t151) * t102) * t96, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:10
	% DurationCPUTime: 1.29s
	% Computational Cost: add. (1304->107), mult. (3345->233), div. (468->12), fcn. (4023->11), ass. (0->103)
	t156 = sin(qJ(2));
	t147 = t156 ^ 2;
	t159 = cos(qJ(2));
	t150 = 0.1e1 / t159 ^ 2;
	t204 = t147 * t150;
	t157 = sin(qJ(1));
	t222 = 0.2e1 * t157;
	t221 = t156 * t204;
	t155 = sin(qJ(6));
	t158 = cos(qJ(6));
	t160 = cos(qJ(1));
	t192 = qJD(6) * t160;
	t197 = qJD(1) * t157;
	t220 = t155 * t197 - t158 * t192;
	t219 = t155 * t192 + t158 * t197;
	t153 = sin(pkin(9));
	t154 = cos(pkin(9));
	t198 = t159 * t160;
	t138 = t153 * t198 - t157 * t154;
	t139 = t157 * t153 + t154 * t198;
	t122 = t138 * t158 + t139 * t155;
	t116 = 0.1e1 / t122;
	t200 = t157 * t156;
	t142 = atan2(-t200, -t159);
	t141 = cos(t142);
	t140 = sin(t142);
	t186 = t140 * t200;
	t126 = -t141 * t159 - t186;
	t123 = 0.1e1 / t126;
	t149 = 0.1e1 / t159;
	t117 = 0.1e1 / t122 ^ 2;
	t124 = 0.1e1 / t126 ^ 2;
	t218 = -0.2e1 * t156;
	t148 = t157 ^ 2;
	t145 = t148 * t204 + 0.1e1;
	t143 = 0.1e1 / t145;
	t217 = t143 - 0.1e1;
	t199 = t157 * t159;
	t136 = -t153 * t199 - t154 * t160;
	t193 = qJD(2) * t160;
	t180 = t156 * t193;
	t129 = t136 * qJD(1) - t153 * t180;
	t137 = t153 * t160 - t154 * t199;
	t130 = t137 * qJD(1) - t154 * t180;
	t105 = t122 * qJD(6) + t129 * t155 - t130 * t158;
	t121 = t138 * t155 - t139 * t158;
	t115 = t121 ^ 2;
	t110 = t115 * t117 + 0.1e1;
	t211 = t117 * t121;
	t106 = -t121 * qJD(6) + t129 * t158 + t130 * t155;
	t118 = t116 * t117;
	t214 = t106 * t118;
	t216 = (t105 * t211 - t115 * t214) / t110 ^ 2;
	t196 = qJD(1) * t160;
	t183 = t156 * t196;
	t194 = qJD(2) * t159;
	t195 = qJD(2) * t157;
	t111 = (-(-t157 * t194 - t183) * t149 + t195 * t204) * t143;
	t206 = t141 * t156;
	t102 = (-t111 * t157 + qJD(2)) * t206 + (-t183 + (t111 - t195) * t159) * t140;
	t215 = t102 * t123 * t124;
	t213 = t111 * t140;
	t212 = t111 * t156;
	t170 = -t153 * t158 - t154 * t155;
	t201 = t156 * t160;
	t134 = t170 * t201;
	t210 = t117 * t134;
	t209 = t124 * t156;
	t202 = t149 * t156;
	t168 = qJD(2) * (t149 * t221 + t202);
	t172 = t147 * t157 * t196;
	t208 = (t148 * t168 + t150 * t172) / t145 ^ 2;
	t177 = 0.1e1 + t204;
	t128 = t177 * t157 * t143;
	t207 = t128 * t157;
	t205 = t147 * t149;
	t152 = t160 ^ 2;
	t203 = t147 * t152;
	t191 = 0.2e1 * t216;
	t114 = t124 * t203 + 0.1e1;
	t190 = 0.2e1 * (-t203 * t215 + (t152 * t156 * t194 - t172) * t124) / t114 ^ 2;
	t189 = 0.2e1 * t215;
	t188 = 0.2e1 * t118 * t121;
	t187 = t124 * t201;
	t185 = t143 * t205;
	t181 = t156 * t195;
	t176 = t156 * t190;
	t175 = t208 * t218;
	t174 = t208 * t222;
	t173 = t157 * t185;
	t171 = t177 * t160;
	t120 = t136 * t158 + t137 * t155;
	t119 = t136 * t155 - t137 * t158;
	t169 = -t153 * t155 + t154 * t158;
	t133 = t169 * t201;
	t132 = -t139 * qJD(1) + t154 * t181;
	t131 = -t138 * qJD(1) + t153 * t181;
	t112 = 0.1e1 / t114;
	t108 = 0.1e1 / t110;
	t107 = (t217 * t156 * t140 - t141 * t173) * t160;
	t104 = -t140 * t199 + t206 + (t140 * t159 - t141 * t200) * t128;
	t103 = -t177 * t174 + (qJD(1) * t171 + t168 * t222) * t143;
	t1 = [t149 * t160 * t175 + (qJD(2) * t171 - t197 * t202) * t143, t103, 0, 0, 0, 0; (t123 * t176 + (-t123 * t194 + (qJD(1) * t107 + t102) * t209) * t112) * t157 + (t124 * t176 * t107 + (-((t111 * t173 + t217 * t194 + t175) * t140 + (t174 * t205 - t212 + (t212 + (t218 - t221) * t195) * t143) * t141) * t187 + (-t124 * t194 + t156 * t189) * t107 + (-t123 + ((-t148 + t152) * t141 * t185 + t217 * t186) * t124) * t156 * qJD(1)) * t112) * t160, (t104 * t209 - t123 * t159) * t160 * t190 + ((-t123 * t197 + (-qJD(2) * t104 - t102) * t160 * t124) * t159 + (-t123 * t193 - (-t103 * t141 * t157 + t140 * t195 + t207 * t213 - t213 + (-qJD(2) * t140 - t141 * t196) * t128) * t187 + (t124 * t197 + t160 * t189) * t104 - ((t103 - t196) * t140 + ((0.1e1 - t207) * qJD(2) + (t128 - t157) * t111) * t141) * t124 * t198) * t156) * t112, 0, 0, 0, 0; (-t116 * t119 + t120 * t211) * t191 + ((t120 * qJD(6) + t131 * t155 - t132 * t158) * t116 + t120 * t106 * t188 + (-t119 * t106 - (-t119 * qJD(6) + t131 * t158 + t132 * t155) * t121 - t120 * t105) * t117) * t108, (-t116 * t133 + t121 * t210) * t191 + (-t105 * t210 + (-t117 * t133 + t134 * t188) * t106 + (t169 * t116 - t170 * t211) * t159 * t193 + ((-t219 * t116 - t220 * t211) * t154 + (t220 * t116 - t219 * t211) * t153) * t156) * t108, 0, 0, 0, -0.2e1 * t216 + 0.2e1 * (t105 * t108 * t117 + (-t108 * t214 - t117 * t216) * t121) * t121;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end