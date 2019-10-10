% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPPR4
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
%   Wie in S6RRPPPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:23
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.79s
	% Computational Cost: add. (774->70), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->75)
	t87 = cos(qJ(1));
	t116 = qJD(1) * t87;
	t85 = sin(qJ(1));
	t138 = 0.2e1 * t85;
	t74 = t85 ^ 2;
	t75 = 0.1e1 / t85;
	t83 = t87 ^ 2;
	t136 = (0.1e1 + 0.1e1 / t74 * t83) * t75 * t116;
	t84 = sin(qJ(2));
	t118 = t85 * t84;
	t86 = cos(qJ(2));
	t65 = atan2(-t118, -t86);
	t63 = sin(t65);
	t108 = t63 * t118;
	t64 = cos(t65);
	t59 = -t64 * t86 - t108;
	t56 = 0.1e1 / t59;
	t79 = 0.1e1 / t86;
	t57 = 0.1e1 / t59 ^ 2;
	t135 = -0.2e1 * t84;
	t73 = t84 ^ 2;
	t80 = 0.1e1 / t86 ^ 2;
	t123 = t73 * t80;
	t70 = t74 * t123 + 0.1e1;
	t66 = 0.1e1 / t70;
	t134 = t66 - 0.1e1;
	t106 = t84 * t116;
	t115 = qJD(2) * t85;
	t125 = t64 * t84;
	t114 = qJD(2) * t86;
	t52 = (-(-t85 * t114 - t106) * t79 + t115 * t123) * t66;
	t48 = (-t52 * t85 + qJD(2)) * t125 + (-t106 + (t52 - t115) * t86) * t63;
	t133 = t48 * t56 * t57;
	t132 = t52 * t63;
	t131 = t52 * t84;
	t130 = t57 * t84;
	t129 = t57 * t87;
	t120 = t79 * t84;
	t72 = t84 * t73;
	t78 = t86 ^ 2;
	t94 = qJD(2) * (t72 * t79 / t78 + t120);
	t99 = t73 * t85 * t116;
	t128 = (t74 * t94 + t80 * t99) / t70 ^ 2;
	t105 = 0.1e1 + t123;
	t62 = t105 * t85 * t66;
	t127 = t62 * t85;
	t126 = t63 * t86;
	t124 = t73 * t79;
	t122 = t73 * t83;
	t76 = 0.1e1 / t85 ^ 2;
	t121 = t76 * t83;
	t119 = t84 * t87;
	t117 = qJD(1) * t85;
	t113 = qJD(2) * t87;
	t55 = t57 * t122 + 0.1e1;
	t98 = t83 * t84 * t114;
	t112 = 0.2e1 * (-t122 * t133 + (t98 - t99) * t57) / t55 ^ 2;
	t111 = 0.2e1 * t133;
	t71 = t78 * t121 + 0.1e1;
	t110 = 0.2e1 * (-t78 * t136 - t76 * t98) / t71 ^ 2;
	t109 = t57 * t119;
	t107 = t66 * t124;
	t104 = 0.1e1 + t121;
	t103 = t84 * t112;
	t102 = t128 * t135;
	t101 = t128 * t138;
	t100 = t85 * t107;
	t97 = t105 * t87;
	t96 = t104 * t84;
	t68 = 0.1e1 / t71;
	t53 = 0.1e1 / t55;
	t51 = (t134 * t84 * t63 - t64 * t100) * t87;
	t50 = -t85 * t126 + t125 + (-t64 * t118 + t126) * t62;
	t49 = -t105 * t101 + (qJD(1) * t97 + t94 * t138) * t66;
	t1 = [t87 * t79 * t102 + (qJD(2) * t97 - t117 * t120) * t66, t49, 0, 0, 0, 0; (t56 * t103 + (-t56 * t114 + (qJD(1) * t51 + t48) * t130) * t53) * t85 + (t57 * t103 * t51 + (-((t52 * t100 + t134 * t114 + t102) * t63 + (t101 * t124 - t131 + (t131 + (-t72 * t80 + t135) * t115) * t66) * t64) * t109 + (t84 * t111 - t57 * t114) * t51 + (-t56 + ((-t74 + t83) * t64 * t107 + t134 * t108) * t57) * t84 * qJD(1)) * t53) * t87, (t50 * t130 - t56 * t86) * t87 * t112 + ((-t56 * t117 + (-qJD(2) * t50 - t48) * t129) * t86 + (-t56 * t113 - (-t49 * t64 * t85 + t63 * t115 + t127 * t132 - t132 + (-qJD(2) * t63 - t116 * t64) * t62) * t109 + (t87 * t111 + t57 * t117) * t50 - ((t49 - t116) * t63 + ((0.1e1 - t127) * qJD(2) + (t62 - t85) * t52) * t64) * t86 * t129) * t84) * t53, 0, 0, 0, 0; t104 * t86 * t110 + (qJD(2) * t96 + 0.2e1 * t86 * t136) * t68, t75 * t110 * t119 + (-t75 * t86 * t113 + qJD(1) * t96) * t68, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:25
	% DurationCPUTime: 0.94s
	% Computational Cost: add. (703->81), mult. (2191->188), div. (456->12), fcn. (2616->9), ass. (0->87)
	t101 = sin(qJ(2));
	t93 = 0.1e1 / t101 ^ 2;
	t103 = cos(qJ(2));
	t97 = t103 ^ 2;
	t151 = t93 * t97;
	t102 = sin(qJ(1));
	t126 = 0.1e1 + t151;
	t95 = t102 ^ 2;
	t91 = t95 * t151 + 0.1e1;
	t89 = 0.1e1 / t91;
	t116 = t126 * t89;
	t75 = t102 * t116;
	t161 = t102 * t75 - 0.1e1;
	t159 = t103 * t151;
	t92 = 0.1e1 / t101;
	t113 = qJD(2) * (-t103 - t159) * t92;
	t104 = cos(qJ(1));
	t138 = qJD(1) * t104;
	t117 = t102 * t97 * t138;
	t160 = (t95 * t113 + t93 * t117) / t91 ^ 2;
	t142 = t102 * t103;
	t88 = atan2(-t142, t101);
	t86 = sin(t88);
	t127 = t86 * t142;
	t87 = cos(t88);
	t70 = t87 * t101 - t127;
	t67 = 0.1e1 / t70;
	t100 = cos(pkin(9));
	t143 = t102 * t100;
	t144 = t101 * t104;
	t99 = sin(pkin(9));
	t83 = t99 * t144 + t143;
	t79 = 0.1e1 / t83;
	t68 = 0.1e1 / t70 ^ 2;
	t80 = 0.1e1 / t83 ^ 2;
	t158 = t89 - 0.1e1;
	t136 = qJD(2) * t102;
	t149 = t101 * t86;
	t63 = ((t101 * t136 - t103 * t138) * t92 + t136 * t151) * t89;
	t58 = (-t63 + t136) * t149 + (-t86 * t138 + (-t102 * t63 + qJD(2)) * t87) * t103;
	t157 = t58 * t67 * t68;
	t98 = t104 ^ 2;
	t150 = t97 * t98;
	t66 = t68 * t150 + 0.1e1;
	t64 = 0.1e1 / t66;
	t155 = t64 * t68;
	t141 = t104 * t100;
	t147 = t102 * t99;
	t114 = t101 * t141 - t147;
	t154 = t80 * t114;
	t153 = t80 * t99;
	t152 = t92 * t97;
	t146 = t104 * t68;
	t145 = qJD(2) * t75;
	t140 = qJD(1) * t102;
	t139 = qJD(1) * t103;
	t137 = qJD(2) * t101;
	t135 = qJD(2) * t103;
	t134 = qJD(2) * t104;
	t133 = 0.2e1 * (-t150 * t157 + (-t101 * t98 * t135 - t117) * t68) / t66 ^ 2;
	t132 = 0.2e1 * t157;
	t78 = t114 ^ 2;
	t73 = t78 * t80 + 0.1e1;
	t122 = t103 * t134;
	t84 = t101 * t143 + t104 * t99;
	t76 = t84 * qJD(1) - t100 * t122;
	t85 = -t101 * t147 + t141;
	t77 = t85 * qJD(1) + t99 * t122;
	t81 = t79 * t80;
	t131 = 0.2e1 * (-t78 * t81 * t77 - t76 * t154) / t73 ^ 2;
	t130 = 0.2e1 * t160;
	t129 = -0.2e1 * t81 * t114;
	t128 = t102 * t89 * t92;
	t125 = t103 * t158;
	t124 = t64 * t137;
	t123 = t102 * t135;
	t121 = t67 * t133;
	t120 = t68 * t133;
	t119 = t103 * t130;
	t118 = t97 * t128;
	t115 = t100 * t79 - t114 * t153;
	t71 = 0.1e1 / t73;
	t112 = t115 * t71;
	t62 = (t87 * t118 + t86 * t125) * t104;
	t60 = -t161 * t87 * t103 + (t102 - t75) * t149;
	t59 = t116 * t138 + 0.2e1 * (t113 * t89 - t126 * t160) * t102;
	t1 = [t128 * t139 + (qJD(2) * t116 + t92 * t119) * t104, t59, 0, 0, 0, 0; (t67 * t124 + (t121 + (qJD(1) * t62 + t58) * t155) * t103) * t102 + (t62 * t68 * t124 + (t62 * t120 + (t62 * t132 + ((t63 * t118 + t158 * t137 + t119) * t86 + (-t63 * t125 + (t130 * t152 + (0.2e1 * t103 + t159) * t89 * qJD(2)) * t102) * t87) * t146) * t64) * t103 + (-t67 + (-(-t95 + t98) * t89 * t87 * t152 + t158 * t127) * t68) * t64 * t139) * t104, (t67 * t64 * t140 + (t121 + (qJD(2) * t60 + t58) * t155) * t104) * t101 + (t60 * t104 * t120 + (-t67 * t134 - ((-t102 * t59 - t138 * t75) * t87 + (t161 * t63 + t136 - t145) * t86) * t103 * t146 + (t104 * t132 + t68 * t140) * t60) * t64 - ((-t59 + t138) * t86 + (-t63 * t75 - qJD(2) + (t63 + t145) * t102) * t87) * t144 * t155) * t103, 0, 0, 0, 0; (-t85 * t154 - t79 * t84) * t131 + ((t114 * qJD(1) + t100 * t123) * t79 + t85 * t77 * t129 + (-t84 * t77 + (-t83 * qJD(1) - t99 * t123) * t114 - t85 * t76) * t80) * t71, t101 * t112 * t134 + (t112 * t140 + (t115 * t131 + (-t76 * t153 + (t100 * t80 + t99 * t129) * t77) * t71) * t104) * t103, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:25
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (926->75), mult. (3228->190), div. (613->15), fcn. (4191->9), ass. (0->83)
	t119 = sin(pkin(9));
	t121 = sin(qJ(2));
	t124 = cos(qJ(1));
	t120 = cos(pkin(9));
	t122 = sin(qJ(1));
	t151 = t122 * t120;
	t106 = t119 * t124 + t121 * t151;
	t123 = cos(qJ(2));
	t150 = t123 * t120;
	t94 = atan2(t106, t150);
	t90 = sin(t94);
	t91 = cos(t94);
	t86 = t106 * t90 + t91 * t150;
	t83 = 0.1e1 / t86;
	t152 = t121 * t124;
	t137 = t120 * t152;
	t153 = t119 * t122;
	t104 = -t137 + t153;
	t115 = 0.1e1 / t123;
	t112 = 0.1e1 / t120;
	t157 = t106 * t112;
	t138 = t115 * t157;
	t103 = t106 ^ 2;
	t113 = 0.1e1 / t120 ^ 2;
	t116 = 0.1e1 / t123 ^ 2;
	t97 = t103 * t113 * t116 + 0.1e1;
	t92 = 0.1e1 / t97;
	t131 = -t90 + (-t91 * t138 + t90) * t92;
	t75 = t131 * t104;
	t169 = 0.2e1 * t75;
	t105 = t119 * t152 + t151;
	t100 = 0.1e1 / t105;
	t101 = 0.1e1 / t105 ^ 2;
	t84 = 0.1e1 / t86 ^ 2;
	t145 = qJD(2) * t124;
	t135 = t123 * t145;
	t88 = t106 * qJD(1) - t120 * t135;
	t165 = t84 * t88;
	t147 = qJD(2) * t121;
	t159 = t123 * t90;
	t161 = t106 * t91;
	t136 = t116 * t147;
	t146 = qJD(2) * t123;
	t149 = qJD(1) * t122;
	t87 = -qJD(1) * t137 + t119 * t149 - t146 * t151;
	t76 = (t106 * t136 - t115 * t87) * t92 * t112;
	t73 = t76 * t161 - t87 * t90 + (-t91 * t147 - t76 * t159) * t120;
	t167 = t73 * t83 * t84;
	t99 = t104 ^ 2;
	t80 = t84 * t99 + 0.1e1;
	t168 = (t104 * t165 - t99 * t167) / t80 ^ 2;
	t114 = t123 ^ 2;
	t117 = t115 / t114;
	t166 = (t103 * t117 * t147 - t106 * t116 * t87) * t113 / t97 ^ 2;
	t107 = t120 * t124 - t121 * t153;
	t133 = t119 * t135;
	t89 = t107 * qJD(1) + t133;
	t164 = t100 * t101 * t89;
	t163 = t104 * t90;
	t162 = t104 * t91;
	t160 = t115 * t92;
	t154 = t116 * t121;
	t82 = (t154 * t157 + t122) * t92;
	t158 = t122 - t82;
	t156 = t107 * t124;
	t118 = t124 ^ 2;
	t155 = t114 * t118;
	t148 = qJD(1) * t124;
	t130 = -t114 * t122 * t148 - t118 * t121 * t146;
	t134 = t155 * t164;
	t139 = t101 * t155;
	t98 = 0.1e1 + t139;
	t144 = 0.2e1 * (t130 * t101 - t134) / t98 ^ 2;
	t143 = -0.2e1 * t166;
	t142 = t83 * t168;
	t141 = t84 * t168;
	t140 = t84 * t163;
	t132 = 0.2e1 * t115 * t166 - t92 * t136;
	t95 = 0.1e1 / t98;
	t78 = 0.1e1 / t80;
	t74 = t82 * t161 + (-t121 * t91 + t158 * t159) * t120;
	t72 = t122 * t143 + t92 * t148 + (-t87 * t92 * t154 + (t143 * t154 + (0.2e1 * t117 * t121 ^ 2 + t115) * t92 * qJD(2)) * t106) * t112;
	t1 = [(t132 * t104 - t88 * t160) * t112, t72, 0, 0, 0, 0; -0.2e1 * t106 * t142 + (-t87 * t83 + (-t106 * t73 - t75 * t88) * t84) * t78 + (t141 * t169 + (t167 * t169 - (t76 * t92 * t138 + t143) * t140 - ((t92 - 0.1e1) * t76 + (t132 * t106 + t87 * t160) * t112) * t84 * t162 - t131 * t165) * t78) * t104, -t74 * t78 * t165 + (-(-t82 * t91 * t87 + (-t76 * t82 * t90 + t72 * t91) * t106) * t84 * t78 + 0.2e1 * (t78 * t167 + t141) * t74) * t104 + (0.2e1 * t124 * t142 * t123 + ((t83 * t145 - (-t158 * qJD(2) + t76) * t140) * t121 + (t83 * t149 + (t124 * t73 - (-t72 + t148) * t163 - (t158 * t76 - qJD(2)) * t162) * t84) * t123) * t78) * t120, 0, 0, 0, 0; (-t100 * t122 - t101 * t156) * t123 * t144 + (-0.2e1 * t123 * t156 * t164 + (-t122 * t147 + t123 * t148) * t100 + (((-t89 - t133) * t122 - t105 * t148) * t123 + (-t121 * t145 - t123 * t149) * t107) * t101) * t95, (-t100 * t152 - t119 * t139) * t144 + (-0.2e1 * t119 * t134 + (-t121 * t149 + t135) * t100 + (0.2e1 * t119 * t130 - t89 * t152) * t101) * t95, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:24
	% EndTime: 2019-10-10 09:23:25
	% DurationCPUTime: 1.24s
	% Computational Cost: add. (1115->105), mult. (3345->233), div. (468->12), fcn. (4023->11), ass. (0->105)
	t158 = sin(qJ(2));
	t149 = 0.1e1 / t158 ^ 2;
	t161 = cos(qJ(2));
	t153 = t161 ^ 2;
	t212 = t149 * t153;
	t229 = t161 * t212;
	t159 = sin(qJ(1));
	t151 = t159 ^ 2;
	t147 = t151 * t212 + 0.1e1;
	t145 = 0.1e1 / t147;
	t148 = 0.1e1 / t158;
	t162 = cos(qJ(1));
	t202 = qJD(1) * t162;
	t187 = t161 * t202;
	t200 = qJD(2) * t159;
	t113 = (-(-t158 * t200 + t187) * t148 + t200 * t212) * t145;
	t228 = t113 - t200;
	t157 = sin(qJ(6));
	t160 = cos(qJ(6));
	t197 = qJD(6) * t162;
	t204 = qJD(1) * t159;
	t227 = -t157 * t204 + t160 * t197;
	t226 = -t157 * t197 - t160 * t204;
	t155 = sin(pkin(9));
	t206 = t162 * t155;
	t156 = cos(pkin(9));
	t209 = t159 * t156;
	t139 = t158 * t206 + t209;
	t205 = t162 * t156;
	t210 = t159 * t155;
	t171 = t158 * t205 - t210;
	t122 = t139 * t160 - t157 * t171;
	t118 = 0.1e1 / t122;
	t208 = t159 * t161;
	t144 = atan2(t208, -t158);
	t143 = cos(t144);
	t142 = sin(t144);
	t191 = t142 * t208;
	t128 = -t143 * t158 + t191;
	t125 = 0.1e1 / t128;
	t119 = 0.1e1 / t122 ^ 2;
	t126 = 0.1e1 / t128 ^ 2;
	t225 = t145 - 0.1e1;
	t140 = t158 * t209 + t206;
	t198 = qJD(2) * t162;
	t185 = t161 * t198;
	t133 = t140 * qJD(1) - t156 * t185;
	t141 = -t158 * t210 + t205;
	t134 = t141 * qJD(1) + t155 * t185;
	t107 = t122 * qJD(6) - t133 * t160 + t134 * t157;
	t175 = -t139 * t157 - t160 * t171;
	t117 = t175 ^ 2;
	t112 = t117 * t119 + 0.1e1;
	t219 = t119 * t175;
	t108 = t175 * qJD(6) + t133 * t157 + t134 * t160;
	t120 = t118 * t119;
	t221 = t108 * t120;
	t224 = (-t107 * t219 - t117 * t221) / t112 ^ 2;
	t154 = t162 ^ 2;
	t211 = t153 * t154;
	t116 = t126 * t211 + 0.1e1;
	t177 = t153 * t159 * t202;
	t199 = qJD(2) * t161;
	t214 = t143 * t161;
	t104 = (t113 * t159 - qJD(2)) * t214 + (t228 * t158 + t187) * t142;
	t222 = t104 * t125 * t126;
	t223 = (-t211 * t222 + (-t154 * t158 * t199 - t177) * t126) / t116 ^ 2;
	t220 = t113 * t161;
	t173 = t155 * t160 - t156 * t157;
	t207 = t161 * t162;
	t136 = t173 * t207;
	t218 = t119 * t136;
	t217 = t126 * t161;
	t216 = t126 * t162;
	t215 = t142 * t159;
	t213 = t148 * t153;
	t203 = qJD(1) * t161;
	t201 = qJD(2) * t158;
	t196 = 0.2e1 * t224;
	t195 = -0.2e1 * t222;
	t194 = -0.2e1 * t120 * t175;
	t170 = qJD(2) * (-t161 - t229) * t148;
	t193 = 0.2e1 * (t149 * t177 + t151 * t170) / t147 ^ 2;
	t192 = t126 * t207;
	t190 = t145 * t213;
	t186 = t159 * t199;
	t182 = 0.1e1 + t212;
	t181 = -0.2e1 * t161 * t223;
	t180 = t159 * t193;
	t179 = t161 * t193;
	t178 = t159 * t190;
	t176 = t182 * t162;
	t174 = t140 * t160 - t141 * t157;
	t124 = t140 * t157 + t141 * t160;
	t172 = t155 * t157 + t156 * t160;
	t135 = t172 * t207;
	t132 = -t139 * qJD(1) - t155 * t186;
	t131 = t171 * qJD(1) + t156 * t186;
	t130 = t182 * t159 * t145;
	t114 = 0.1e1 / t116;
	t111 = (-t225 * t161 * t142 - t143 * t178) * t162;
	t109 = 0.1e1 / t112;
	t106 = -t158 * t215 - t214 + (t142 * t158 + t143 * t208) * t130;
	t105 = -t182 * t180 + (qJD(1) * t176 + 0.2e1 * t159 * t170) * t145;
	t1 = [t162 * t148 * t179 + (t148 * t159 * t203 + qJD(2) * t176) * t145, t105, 0, 0, 0, 0; (t125 * t181 + (-t125 * t201 + (-qJD(1) * t111 - t104) * t217) * t114) * t159 + (t126 * t181 * t111 + (((t113 * t178 + t225 * t201 + t179) * t142 + (t180 * t213 + t220 + (-t220 + (0.2e1 * t161 + t229) * t200) * t145) * t143) * t192 + (-t126 * t201 + t161 * t195) * t111 + (t125 + ((t151 - t154) * t143 * t190 + t225 * t191) * t126) * t203) * t114) * t162, 0.2e1 * (-t106 * t217 - t125 * t158) * t162 * t223 + ((-t125 * t204 + (-qJD(2) * t106 - t104) * t216) * t158 + (t125 * t198 + (t105 * t143 * t159 + t228 * t142 + (qJD(2) * t142 - t113 * t215 + t143 * t202) * t130) * t192 + (-t126 * t204 + t162 * t195) * t106 + ((t105 - t202) * t142 + ((-t130 * t159 + 0.1e1) * qJD(2) + (t130 - t159) * t113) * t143) * t158 * t216) * t161) * t114, 0, 0, 0, 0; (t118 * t174 - t124 * t219) * t196 + ((t124 * qJD(6) - t131 * t160 + t132 * t157) * t118 + t124 * t108 * t194 + (t174 * t108 + (t174 * qJD(6) + t131 * t157 + t132 * t160) * t175 - t124 * t107) * t119) * t109, (-t118 * t135 - t175 * t218) * t196 + (-t107 * t218 + (-t135 * t119 + t136 * t194) * t108 + (-t172 * t118 - t173 * t219) * t158 * t198 + ((t226 * t118 - t227 * t219) * t156 + (t227 * t118 + t226 * t219) * t155) * t161) * t109, 0, 0, 0, -0.2e1 * t224 - 0.2e1 * (t107 * t119 * t109 - (-t109 * t221 - t119 * t224) * t175) * t175;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end