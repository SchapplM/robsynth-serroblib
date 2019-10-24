% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPP2
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5PRRPP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:29
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRPP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:09
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:09
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:09
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:09
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (429->41), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->59)
	t97 = sin(qJ(2));
	t90 = t97 ^ 2;
	t99 = cos(qJ(2));
	t123 = t90 / t99 ^ 2;
	t111 = 0.1e1 + t123;
	t94 = sin(pkin(7));
	t135 = t111 * t94;
	t134 = qJD(2) * t97;
	t95 = cos(pkin(7));
	t120 = t95 * t99;
	t96 = sin(qJ(3));
	t98 = cos(qJ(3));
	t77 = t98 * t120 + t94 * t96;
	t121 = t94 * t97;
	t80 = atan2(-t121, -t99);
	t78 = sin(t80);
	t79 = cos(t80);
	t63 = -t78 * t121 - t79 * t99;
	t60 = 0.1e1 / t63;
	t73 = 0.1e1 / t77;
	t124 = t78 * t99;
	t108 = -t94 * t124 + t79 * t97;
	t115 = t79 * t121;
	t109 = -t115 + t124;
	t87 = t94 ^ 2;
	t83 = t87 * t123 + 0.1e1;
	t81 = 0.1e1 / t83;
	t67 = t81 * t135;
	t54 = t109 * t67 + t108;
	t132 = 0.2e1 * t54;
	t61 = 0.1e1 / t63 ^ 2;
	t74 = 0.1e1 / t77 ^ 2;
	t118 = qJD(2) * t99;
	t127 = t61 * t97;
	t64 = qJD(2) * t67;
	t53 = t108 * qJD(2) + t109 * t64;
	t130 = t53 * t60 * t61;
	t88 = t95 ^ 2;
	t59 = t88 * t90 * t61 + 0.1e1;
	t131 = (t118 * t127 - t90 * t130) * t88 / t59 ^ 2;
	t76 = t96 * t120 - t94 * t98;
	t125 = t74 * t76;
	t113 = t95 * t134;
	t117 = qJD(3) * t76;
	t70 = -t98 * t113 - t117;
	t126 = t70 * t73 * t74;
	t72 = t76 ^ 2;
	t68 = t72 * t74 + 0.1e1;
	t69 = -t77 * qJD(3) + t96 * t113;
	t129 = (-t69 * t125 - t72 * t126) / t68 ^ 2;
	t57 = 0.1e1 / t59;
	t128 = t57 * t61;
	t119 = t67 - t94;
	t116 = -0.2e1 * t129;
	t112 = t67 * t94 - 0.1e1;
	t110 = t98 * t125 - t73 * t96;
	t65 = 0.1e1 / t68;
	t56 = 0.2e1 * (t81 - t111 / t83 ^ 2 * t87) / t99 * t134 * t135;
	t1 = [0, t56, 0, 0, 0; 0, ((-0.2e1 * t60 * t131 + (-qJD(2) * t54 - t53) * t128) * t99 + (t61 * t131 * t132 + (t56 * t61 * t115 + t130 * t132 - qJD(2) * t60 - (-t119 * qJD(2) + t112 * t64) * t78 * t127) * t57 - (t56 * t78 + (-t112 * qJD(2) + t119 * t64) * t79) * t99 * t128) * t97) * t95, 0, 0, 0; 0, (t110 * t97 * t116 + (t110 * t118 + ((-qJD(3) * t73 - 0.2e1 * t76 * t126) * t98 + (-t69 * t98 + (t70 - t117) * t96) * t74) * t97) * t65) * t95, t116 + 0.2e1 * (-t65 * t69 * t74 + (-t65 * t126 - t74 * t129) * t76) * t76, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:09
	% EndTime: 2019-10-24 10:29:09
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (614->45), mult. (1166->104), div. (229->12), fcn. (1315->9), ass. (0->57)
	t105 = sin(pkin(7));
	t107 = sin(qJ(2));
	t101 = t107 ^ 2;
	t108 = cos(qJ(2));
	t126 = t101 / t108 ^ 2;
	t118 = 0.1e1 + t126;
	t141 = t105 * t118;
	t140 = qJD(2) * t107;
	t106 = cos(pkin(7));
	t124 = t106 * t108;
	t99 = qJ(3) + pkin(8);
	t95 = sin(t99);
	t96 = cos(t99);
	t85 = t105 * t95 + t96 * t124;
	t125 = t105 * t107;
	t89 = atan2(-t125, -t108);
	t86 = sin(t89);
	t87 = cos(t89);
	t74 = -t87 * t108 - t86 * t125;
	t71 = 0.1e1 / t74;
	t81 = 0.1e1 / t85;
	t128 = t87 * t107;
	t129 = t108 * t86;
	t97 = t105 ^ 2;
	t92 = t97 * t126 + 0.1e1;
	t90 = 0.1e1 / t92;
	t78 = t90 * t141;
	t63 = t78 * t129 + t128 + (-t78 * t128 - t129) * t105;
	t138 = 0.2e1 * t63;
	t72 = 0.1e1 / t74 ^ 2;
	t82 = 0.1e1 / t85 ^ 2;
	t123 = qJD(2) * t108;
	t130 = t107 * t72;
	t121 = t87 * t125;
	t75 = qJD(2) * t78;
	t62 = (-t121 + t129) * t75 + (-t105 * t129 + t128) * qJD(2);
	t135 = t62 * t71 * t72;
	t98 = t106 ^ 2;
	t67 = t98 * t101 * t72 + 0.1e1;
	t137 = (-t101 * t135 + t123 * t130) * t98 / t67 ^ 2;
	t84 = -t105 * t96 + t95 * t124;
	t132 = t82 * t84;
	t119 = t106 * t140;
	t127 = qJD(3) * t84;
	t77 = -t96 * t119 - t127;
	t133 = t77 * t81 * t82;
	t80 = t84 ^ 2;
	t70 = t80 * t82 + 0.1e1;
	t76 = -t85 * qJD(3) + t95 * t119;
	t136 = (-t76 * t132 - t80 * t133) / t70 ^ 2;
	t65 = 0.1e1 / t67;
	t134 = t65 * t72;
	t122 = -0.2e1 * t136;
	t117 = t96 * t132 - t81 * t95;
	t68 = 0.1e1 / t70;
	t64 = 0.2e1 * (t90 - t118 / t92 ^ 2 * t97) / t108 * t140 * t141;
	t1 = [0, t64, 0, 0, 0; 0, ((-0.2e1 * t71 * t137 + (-qJD(2) * t63 - t62) * t134) * t108 + (t72 * t137 * t138 + (t64 * t72 * t121 + t135 * t138 - qJD(2) * t71 - ((t105 * t78 - 0.1e1) * t75 + (t105 - t78) * qJD(2)) * t86 * t130) * t65 - (t64 * t86 + (qJD(2) + (-0.2e1 * t105 + t78) * t75) * t87) * t108 * t134) * t107) * t106, 0, 0, 0; 0, (t117 * t68 * t123 + (t117 * t122 + ((-qJD(3) * t81 - 0.2e1 * t84 * t133) * t96 + (-t76 * t96 + (t77 - t127) * t95) * t82) * t68) * t107) * t106, t122 + 0.2e1 * (-t68 * t76 * t82 + (-t68 * t133 - t82 * t136) * t84) * t84, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:09
	% EndTime: 2019-10-24 10:29:10
	% DurationCPUTime: 0.70s
	% Computational Cost: add. (3088->84), mult. (3936->201), div. (821->15), fcn. (5024->9), ass. (0->85)
	t132 = qJ(3) + pkin(8);
	t129 = sin(t132);
	t130 = cos(t132);
	t138 = cos(pkin(7));
	t137 = sin(pkin(7));
	t140 = cos(qJ(2));
	t167 = t137 * t140;
	t116 = t129 * t167 + t138 * t130;
	t139 = sin(qJ(2));
	t165 = t139 * t129;
	t106 = atan2(-t116, t165);
	t102 = sin(t106);
	t103 = cos(t106);
	t111 = t116 ^ 2;
	t127 = 0.1e1 / t129 ^ 2;
	t135 = 0.1e1 / t139 ^ 2;
	t169 = t127 * t135;
	t107 = t111 * t169 + 0.1e1;
	t104 = 0.1e1 / t107;
	t126 = 0.1e1 / t129;
	t157 = t126 * t135 * t140;
	t150 = t116 * t157 + t137;
	t92 = t150 * t104;
	t176 = t137 - t92;
	t185 = t102 * t139 * t176 + t103 * t140;
	t166 = t138 * t140;
	t120 = t137 * t129 + t130 * t166;
	t114 = 0.1e1 / t120 ^ 2;
	t131 = t138 ^ 2;
	t133 = t139 ^ 2;
	t168 = t131 * t133;
	t110 = t114 * t168 + 0.1e1;
	t163 = qJD(2) * t140;
	t152 = t139 * t163;
	t119 = t129 * t166 - t137 * t130;
	t164 = qJD(2) * t139;
	t153 = t138 * t164;
	t101 = -qJD(3) * t119 - t130 * t153;
	t113 = 0.1e1 / t120;
	t175 = t101 * t113 * t114;
	t184 = 0.1e1 / t110 ^ 2 * (t114 * t152 - t133 * t175) * t131;
	t134 = 0.1e1 / t139;
	t118 = -t138 * t129 + t130 * t167;
	t170 = t127 * t130;
	t149 = t116 * t170 - t118 * t126;
	t183 = t134 * t149;
	t128 = t126 * t127;
	t136 = t134 / t133;
	t158 = t116 * t169;
	t161 = qJD(3) * t130;
	t151 = t140 * t161;
	t154 = t137 * t164;
	t98 = -t137 * t151 + (qJD(3) * t138 + t154) * t129;
	t182 = -0.2e1 / t107 ^ 2 * (-t98 * t158 + (-t127 * t136 * t163 - t128 * t135 * t161) * t111);
	t174 = t102 * t116;
	t97 = t165 * t103 - t174;
	t93 = 0.1e1 / t97;
	t94 = 0.1e1 / t97 ^ 2;
	t181 = 0.2e1 * t119;
	t155 = t129 * t163;
	t148 = t139 * t161 + t155;
	t86 = (t98 * t134 * t126 + t148 * t158) * t104;
	t178 = t116 * t86;
	t83 = (-t86 * t165 + t98) * t102 + (t148 - t178) * t103;
	t95 = t93 * t94;
	t180 = t83 * t95;
	t162 = qJD(3) * t129;
	t100 = t129 * t153 - t137 * t162 - t138 * t151;
	t179 = t100 * t94;
	t177 = t119 * t94;
	t172 = t103 * t116;
	t112 = t119 ^ 2;
	t91 = t112 * t94 + 0.1e1;
	t160 = 0.2e1 * (-t100 * t177 - t112 * t180) / t91 ^ 2;
	t159 = t139 * t181;
	t156 = t130 * t168;
	t108 = 0.1e1 / t110;
	t99 = t116 * qJD(3) + t130 * t154;
	t89 = 0.1e1 / t91;
	t88 = t104 * t183;
	t85 = t185 * t129 - t92 * t172;
	t84 = (-t116 * t88 + t130 * t139) * t103 + (-t88 * t165 - t118) * t102;
	t82 = t150 * t182 + (-t98 * t157 + (-t151 * t169 + (-0.2e1 * t136 * t140 ^ 2 - t134) * t126 * qJD(2)) * t116) * t104;
	t80 = t182 * t183 + (-t149 * t135 * t163 + (-t98 * t170 + t126 * t99 + (t118 * t170 + (-0.2e1 * t128 * t130 ^ 2 - t126) * t116) * qJD(3)) * t134) * t104;
	t1 = [0, t82, t80, 0, 0; 0, t85 * t160 * t177 + (t85 * t179 + (-(-t172 * t82 + (t103 * t98 + t174 * t86) * t92) * t94 + 0.2e1 * t85 * t180) * t119 + (-t138 * t139 * t93 - t185 * t177) * t161) * t89 + (-((t176 * t86 - qJD(2)) * t139 * t103 + (-t139 * t82 + (qJD(2) * t176 - t86) * t140) * t102) * t89 * t177 + (-t93 * t89 * t163 + (t83 * t89 * t94 + t93 * t160) * t139) * t138) * t129, (-t120 * t93 + t177 * t84) * t160 + (t84 * t179 + t101 * t93 + (t181 * t84 * t95 - t120 * t94) * t83 + (-(t130 * t163 - t116 * t80 - t118 * t86 + t88 * t98 + (-t86 * t88 - qJD(3)) * t165) * t103 - (t99 + (-t155 + t178) * t88 + (-t129 * t80 + (-qJD(3) * t88 - t86) * t130) * t139) * t102) * t177) * t89, 0, 0; 0, 0.2e1 * (t113 * t166 + t114 * t156) * t184 + (0.2e1 * t156 * t175 + t113 * t153 + (t101 * t166 + (-0.2e1 * t130 * t152 + t133 * t162) * t131) * t114) * t108, (t114 * t159 * t184 + (t159 * t175 + (t100 * t139 - t119 * t163) * t114) * t108) * t138, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end