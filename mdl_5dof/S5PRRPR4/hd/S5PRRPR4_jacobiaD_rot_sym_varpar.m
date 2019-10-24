% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPR4
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
%   Wie in S5PRRPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:31
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5PRRPR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:00
	% EndTime: 2019-10-24 10:31:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:00
	% EndTime: 2019-10-24 10:31:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:00
	% EndTime: 2019-10-24 10:31:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:00
	% EndTime: 2019-10-24 10:31:00
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (429->41), mult. (1166->102), div. (229->12), fcn. (1315->9), ass. (0->59)
	t97 = sin(qJ(2));
	t90 = t97 ^ 2;
	t99 = cos(qJ(2));
	t123 = t90 / t99 ^ 2;
	t111 = 0.1e1 + t123;
	t94 = sin(pkin(8));
	t135 = t111 * t94;
	t134 = qJD(2) * t97;
	t95 = cos(pkin(8));
	t120 = t95 * t99;
	t96 = sin(qJ(3));
	t98 = cos(qJ(3));
	t77 = t98 * t120 + t94 * t96;
	t121 = t94 * t97;
	t80 = atan2(-t121, -t99);
	t78 = sin(t80);
	t79 = cos(t80);
	t63 = -t121 * t78 - t79 * t99;
	t60 = 0.1e1 / t63;
	t73 = 0.1e1 / t77;
	t124 = t78 * t99;
	t108 = -t124 * t94 + t79 * t97;
	t115 = t79 * t121;
	t109 = -t115 + t124;
	t87 = t94 ^ 2;
	t83 = t123 * t87 + 0.1e1;
	t81 = 0.1e1 / t83;
	t67 = t81 * t135;
	t54 = t109 * t67 + t108;
	t132 = 0.2e1 * t54;
	t61 = 0.1e1 / t63 ^ 2;
	t74 = 0.1e1 / t77 ^ 2;
	t118 = qJD(2) * t99;
	t127 = t61 * t97;
	t64 = qJD(2) * t67;
	t53 = qJD(2) * t108 + t109 * t64;
	t130 = t53 * t60 * t61;
	t88 = t95 ^ 2;
	t59 = t61 * t88 * t90 + 0.1e1;
	t131 = (t118 * t127 - t130 * t90) * t88 / t59 ^ 2;
	t76 = t120 * t96 - t94 * t98;
	t125 = t74 * t76;
	t113 = t95 * t134;
	t117 = qJD(3) * t76;
	t70 = -t113 * t98 - t117;
	t126 = t70 * t73 * t74;
	t72 = t76 ^ 2;
	t68 = t72 * t74 + 0.1e1;
	t69 = -qJD(3) * t77 + t96 * t113;
	t129 = (-t125 * t69 - t126 * t72) / t68 ^ 2;
	t57 = 0.1e1 / t59;
	t128 = t57 * t61;
	t119 = t67 - t94;
	t116 = -0.2e1 * t129;
	t112 = t67 * t94 - 0.1e1;
	t110 = t125 * t98 - t73 * t96;
	t65 = 0.1e1 / t68;
	t56 = 0.2e1 * (t81 - t111 / t83 ^ 2 * t87) / t99 * t134 * t135;
	t1 = [0, t56, 0, 0, 0; 0, ((-0.2e1 * t60 * t131 + (-qJD(2) * t54 - t53) * t128) * t99 + (t61 * t131 * t132 + (t56 * t61 * t115 + t130 * t132 - qJD(2) * t60 - (-qJD(2) * t119 + t112 * t64) * t78 * t127) * t57 - (t56 * t78 + (-qJD(2) * t112 + t119 * t64) * t79) * t99 * t128) * t97) * t95, 0, 0, 0; 0, (t110 * t97 * t116 + (t110 * t118 + ((-qJD(3) * t73 - 0.2e1 * t126 * t76) * t98 + (-t69 * t98 + (t70 - t117) * t96) * t74) * t97) * t65) * t95, t116 + 0.2e1 * (-t65 * t69 * t74 + (-t126 * t65 - t129 * t74) * t76) * t76, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:00
	% EndTime: 2019-10-24 10:31:00
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (614->45), mult. (1166->104), div. (229->12), fcn. (1315->9), ass. (0->57)
	t105 = sin(pkin(8));
	t107 = sin(qJ(2));
	t101 = t107 ^ 2;
	t108 = cos(qJ(2));
	t126 = t101 / t108 ^ 2;
	t118 = 0.1e1 + t126;
	t141 = t105 * t118;
	t140 = qJD(2) * t107;
	t106 = cos(pkin(8));
	t124 = t106 * t108;
	t99 = qJ(3) + pkin(9);
	t95 = sin(t99);
	t96 = cos(t99);
	t85 = t105 * t95 + t96 * t124;
	t125 = t105 * t107;
	t89 = atan2(-t125, -t108);
	t86 = sin(t89);
	t87 = cos(t89);
	t74 = -t108 * t87 - t125 * t86;
	t71 = 0.1e1 / t74;
	t81 = 0.1e1 / t85;
	t128 = t87 * t107;
	t129 = t108 * t86;
	t97 = t105 ^ 2;
	t92 = t126 * t97 + 0.1e1;
	t90 = 0.1e1 / t92;
	t78 = t90 * t141;
	t63 = t78 * t129 + t128 + (-t128 * t78 - t129) * t105;
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
	t67 = t101 * t72 * t98 + 0.1e1;
	t137 = (-t101 * t135 + t123 * t130) * t98 / t67 ^ 2;
	t84 = -t105 * t96 + t124 * t95;
	t132 = t82 * t84;
	t119 = t106 * t140;
	t127 = qJD(3) * t84;
	t77 = -t119 * t96 - t127;
	t133 = t77 * t81 * t82;
	t80 = t84 ^ 2;
	t70 = t80 * t82 + 0.1e1;
	t76 = -qJD(3) * t85 + t95 * t119;
	t136 = (-t132 * t76 - t133 * t80) / t70 ^ 2;
	t65 = 0.1e1 / t67;
	t134 = t65 * t72;
	t122 = -0.2e1 * t136;
	t117 = t132 * t96 - t81 * t95;
	t68 = 0.1e1 / t70;
	t64 = 0.2e1 * (t90 - t118 / t92 ^ 2 * t97) / t108 * t140 * t141;
	t1 = [0, t64, 0, 0, 0; 0, ((-0.2e1 * t71 * t137 + (-qJD(2) * t63 - t62) * t134) * t108 + (t72 * t137 * t138 + (t64 * t72 * t121 + t135 * t138 - qJD(2) * t71 - ((t105 * t78 - 0.1e1) * t75 + (t105 - t78) * qJD(2)) * t86 * t130) * t65 - (t64 * t86 + (qJD(2) + (-0.2e1 * t105 + t78) * t75) * t87) * t108 * t134) * t107) * t106, 0, 0, 0; 0, (t117 * t68 * t123 + (t117 * t122 + ((-qJD(3) * t81 - 0.2e1 * t133 * t84) * t96 + (-t76 * t96 + (t77 - t127) * t95) * t82) * t68) * t107) * t106, t122 + 0.2e1 * (-t68 * t76 * t82 + (-t133 * t68 - t136 * t82) * t84) * t84, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:00
	% EndTime: 2019-10-24 10:31:00
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (1050->45), mult. (1346->105), div. (247->12), fcn. (1511->9), ass. (0->61)
	t127 = sin(pkin(8));
	t129 = sin(qJ(2));
	t123 = t129 ^ 2;
	t130 = cos(qJ(2));
	t151 = t123 / t130 ^ 2;
	t140 = 0.1e1 + t151;
	t165 = t127 * t140;
	t164 = qJD(2) * t129;
	t118 = qJ(3) + pkin(9) + qJ(5);
	t116 = sin(t118);
	t117 = cos(t118);
	t128 = cos(pkin(8));
	t148 = t128 * t130;
	t105 = t127 * t116 + t117 * t148;
	t149 = t127 * t129;
	t111 = atan2(-t149, -t130);
	t108 = sin(t111);
	t109 = cos(t111);
	t97 = -t108 * t149 - t109 * t130;
	t94 = 0.1e1 / t97;
	t119 = t127 ^ 2;
	t114 = t119 * t151 + 0.1e1;
	t112 = 0.1e1 / t114;
	t99 = t112 * t165;
	t142 = t127 * t99 - 0.1e1;
	t152 = t109 * t129;
	t153 = t108 * t130;
	t155 = t127 - t99;
	t84 = -t142 * t152 - t153 * t155;
	t162 = 0.2e1 * t84;
	t101 = 0.1e1 / t105;
	t102 = 0.1e1 / t105 ^ 2;
	t95 = 0.1e1 / t97 ^ 2;
	t120 = t128 ^ 2;
	t147 = qJD(2) * t130;
	t156 = t129 * t95;
	t145 = t109 * t149;
	t98 = qJD(2) * t99;
	t83 = (-t145 + t153) * t98 + (-t127 * t153 + t152) * qJD(2);
	t159 = t83 * t94 * t95;
	t91 = t120 * t123 * t95 + 0.1e1;
	t161 = (-t123 * t159 + t147 * t156) * t120 / t91 ^ 2;
	t144 = t116 * t148;
	t104 = -t117 * t127 + t144;
	t100 = t104 ^ 2;
	t154 = t102 * t104;
	t121 = qJD(3) + qJD(5);
	t141 = t128 * t164;
	t93 = -t121 * t144 + (t121 * t127 - t141) * t117;
	t157 = t101 * t102 * t93;
	t88 = t100 * t102 + 0.1e1;
	t92 = -t105 * t121 + t116 * t141;
	t160 = (-t100 * t157 - t154 * t92) / t88 ^ 2;
	t89 = 0.1e1 / t91;
	t158 = t89 * t95;
	t146 = -0.2e1 * t160;
	t139 = -t101 * t116 + t117 * t154;
	t86 = 0.1e1 / t88;
	t85 = 0.2e1 * (t112 - t140 / t114 ^ 2 * t119) / t130 * t164 * t165;
	t80 = t146 + 0.2e1 * (-t102 * t86 * t92 + (-t102 * t160 - t157 * t86) * t104) * t104;
	t1 = [0, t85, 0, 0, 0; 0, ((-0.2e1 * t94 * t161 + (-qJD(2) * t84 - t83) * t158) * t130 + (t95 * t161 * t162 + (t85 * t95 * t145 + t159 * t162 - qJD(2) * t94 - (qJD(2) * t155 + t142 * t98) * t108 * t156) * t89 - (t108 * t85 + (qJD(2) + (-0.2e1 * t127 + t99) * t98) * t109) * t130 * t158) * t129) * t128, 0, 0, 0; 0, (t139 * t86 * t147 + (t139 * t146 + ((-t101 * t121 - 0.2e1 * t104 * t157) * t117 + (-t117 * t92 + (-t104 * t121 + t93) * t116) * t102) * t86) * t129) * t128, t80, 0, t80;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end