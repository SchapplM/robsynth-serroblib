% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPP3
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
%   Wie in S5RPRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRPP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_jacobiaD_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:13:29
	% EndTime: 2019-12-31 18:13:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:13:29
	% EndTime: 2019-12-31 18:13:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:13:29
	% EndTime: 2019-12-31 18:13:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:13:29
	% EndTime: 2019-12-31 18:13:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:13:29
	% EndTime: 2019-12-31 18:13:30
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (1883->70), mult. (1839->155), div. (436->14), fcn. (2165->7), ass. (0->75)
	t92 = sin(qJ(1));
	t142 = 0.2e1 * t92;
	t87 = 0.1e1 / t92;
	t93 = cos(qJ(1));
	t126 = t87 * t93;
	t120 = qJD(3) * t92;
	t122 = qJD(1) * t93;
	t85 = pkin(7) + qJ(3);
	t83 = sin(t85);
	t112 = t83 * t122;
	t78 = t83 ^ 2;
	t84 = cos(t85);
	t81 = 0.1e1 / t84 ^ 2;
	t129 = t78 * t81;
	t86 = t92 ^ 2;
	t73 = t86 * t129 + 0.1e1;
	t71 = 0.1e1 / t73;
	t80 = 0.1e1 / t84;
	t57 = (-(-t84 * t120 - t112) * t80 + t120 * t129) * t71;
	t141 = t57 - t120;
	t91 = t93 ^ 2;
	t140 = qJD(1) * (0.1e1 / t86 * t91 + 0.1e1) * t126;
	t124 = t92 * t83;
	t70 = atan2(-t124, -t84);
	t68 = sin(t70);
	t114 = t68 * t124;
	t69 = cos(t70);
	t64 = -t69 * t84 - t114;
	t61 = 0.1e1 / t64;
	t62 = 0.1e1 / t64 ^ 2;
	t139 = -0.2e1 * t83;
	t138 = t71 - 0.1e1;
	t131 = t69 * t83;
	t53 = (-t57 * t92 + qJD(3)) * t131 + (t141 * t84 - t112) * t68;
	t137 = t53 * t61 * t62;
	t136 = t57 * t83;
	t135 = t62 * t83;
	t134 = t62 * t93;
	t127 = t80 * t83;
	t77 = t83 * t78;
	t79 = t84 ^ 2;
	t100 = qJD(3) * (t77 * t80 / t79 + t127);
	t105 = t78 * t92 * t122;
	t133 = (t86 * t100 + t81 * t105) / t73 ^ 2;
	t132 = t68 * t92;
	t130 = t78 * t80;
	t128 = t78 * t91;
	t88 = 0.1e1 / t92 ^ 2;
	t125 = t88 * t91;
	t123 = qJD(1) * t92;
	t121 = qJD(3) * t84;
	t119 = qJD(3) * t93;
	t104 = t83 * t91 * t121;
	t60 = t62 * t128 + 0.1e1;
	t118 = 0.2e1 * (-t128 * t137 + (t104 - t105) * t62) / t60 ^ 2;
	t117 = 0.2e1 * t137;
	t76 = t79 * t125 + 0.1e1;
	t116 = 0.2e1 * (-t88 * t104 - t79 * t140) / t76 ^ 2;
	t115 = t83 * t134;
	t113 = t71 * t130;
	t111 = 0.1e1 + t129;
	t110 = 0.1e1 + t125;
	t109 = t83 * t118;
	t108 = t133 * t139;
	t107 = t133 * t142;
	t106 = t92 * t113;
	t103 = t111 * t93;
	t102 = t110 * t83;
	t74 = 0.1e1 / t76;
	t66 = t111 * t92 * t71;
	t58 = 0.1e1 / t60;
	t56 = (t138 * t83 * t68 - t69 * t106) * t93;
	t55 = -t84 * t132 + t131 + (-t69 * t124 + t68 * t84) * t66;
	t54 = -t111 * t107 + (qJD(1) * t103 + t100 * t142) * t71;
	t1 = [t93 * t80 * t108 + (qJD(3) * t103 - t123 * t127) * t71, 0, t54, 0, 0; (t61 * t109 + (-t61 * t121 + (qJD(1) * t56 + t53) * t135) * t58) * t92 + (t62 * t109 * t56 + (-((t57 * t106 + t138 * t121 + t108) * t68 + (t107 * t130 - t136 + (t136 + (-t77 * t81 + t139) * t120) * t71) * t69) * t115 + (t83 * t117 - t62 * t121) * t56 + (-t61 + ((-t86 + t91) * t69 * t113 + t138 * t114) * t62) * t83 * qJD(1)) * t58) * t93, 0, (t55 * t135 - t61 * t84) * t93 * t118 + ((-t61 * t123 + (-qJD(3) * t55 - t53) * t134) * t84 + (-t61 * t119 - (-t54 * t69 * t92 - t141 * t68 + (-qJD(3) * t68 - t122 * t69 + t132 * t57) * t66) * t115 + (t93 * t117 + t62 * t123) * t55 - ((t54 - t122) * t68 + ((-t66 * t92 + 0.1e1) * qJD(3) + (t66 - t92) * t57) * t69) * t84 * t134) * t83) * t58, 0, 0; t110 * t84 * t116 + (qJD(3) * t102 + 0.2e1 * t84 * t140) * t74, 0, t83 * t116 * t126 + (-t84 * t87 * t119 + qJD(1) * t102) * t74, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:13:29
	% EndTime: 2019-12-31 18:13:30
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (1694->69), mult. (1839->161), div. (436->14), fcn. (2165->7), ass. (0->74)
	t92 = sin(qJ(1));
	t87 = 0.1e1 / t92;
	t93 = cos(qJ(1));
	t125 = t87 * t93;
	t86 = t92 ^ 2;
	t91 = t93 ^ 2;
	t101 = qJD(1) * (-0.1e1 / t86 * t91 - 0.1e1) * t125;
	t85 = pkin(7) + qJ(3);
	t83 = sin(t85);
	t120 = qJD(3) * t83;
	t84 = cos(t85);
	t104 = t84 * t91 * t120;
	t88 = 0.1e1 / t92 ^ 2;
	t124 = t88 * t91;
	t77 = t83 ^ 2;
	t76 = t77 * t124 + 0.1e1;
	t140 = (t77 * t101 + t88 * t104) / t76 ^ 2;
	t79 = 0.1e1 / t83 ^ 2;
	t82 = t84 ^ 2;
	t128 = t79 * t82;
	t109 = 0.1e1 + t128;
	t139 = t109 * t92;
	t119 = qJD(3) * t92;
	t121 = qJD(1) * t93;
	t110 = t84 * t121;
	t73 = t86 * t128 + 0.1e1;
	t71 = 0.1e1 / t73;
	t78 = 0.1e1 / t83;
	t57 = ((t83 * t119 - t110) * t78 + t119 * t128) * t71;
	t138 = -t57 + t119;
	t123 = t92 * t84;
	t70 = atan2(-t123, t83);
	t68 = sin(t70);
	t112 = t68 * t123;
	t69 = cos(t70);
	t64 = t69 * t83 - t112;
	t61 = 0.1e1 / t64;
	t62 = 0.1e1 / t64 ^ 2;
	t137 = t71 - 0.1e1;
	t130 = t69 * t84;
	t53 = (-t57 * t92 + qJD(3)) * t130 + (t138 * t83 - t110) * t68;
	t136 = t53 * t61 * t62;
	t135 = t57 * t84;
	t134 = t62 * t84;
	t133 = t62 * t93;
	t129 = t78 * t84;
	t81 = t84 * t82;
	t100 = qJD(3) * (-t78 / t77 * t81 - t129);
	t126 = t82 * t92;
	t105 = t121 * t126;
	t132 = (t86 * t100 + t79 * t105) / t73 ^ 2;
	t131 = t68 * t92;
	t127 = t82 * t91;
	t122 = qJD(1) * t92;
	t118 = qJD(3) * t93;
	t60 = t62 * t127 + 0.1e1;
	t117 = 0.2e1 * (-t127 * t136 + (-t104 - t105) * t62) / t60 ^ 2;
	t116 = 0.2e1 * t136;
	t115 = -0.2e1 * t132;
	t114 = t84 * t133;
	t113 = t84 * t132;
	t111 = t78 * t126;
	t108 = 0.1e1 + t124;
	t107 = t84 * t117;
	t106 = t71 * t111;
	t103 = t109 * t93;
	t102 = t108 * t84;
	t74 = 0.1e1 / t76;
	t66 = t71 * t139;
	t58 = 0.1e1 / t60;
	t56 = (t137 * t84 * t68 + t69 * t106) * t93;
	t55 = t83 * t131 + t130 + (-t69 * t123 - t68 * t83) * t66;
	t54 = t115 * t139 + (qJD(1) * t103 + 0.2e1 * t100 * t92) * t71;
	t1 = [0.2e1 * t93 * t78 * t113 + (qJD(3) * t103 + t122 * t129) * t71, 0, t54, 0, 0; (t61 * t107 + (t61 * t120 + (qJD(1) * t56 + t53) * t134) * t58) * t92 + (t62 * t107 * t56 + (-((-t57 * t106 - t137 * t120 - 0.2e1 * t113) * t68 + (t111 * t115 - t135 + (t135 + (-t79 * t81 - 0.2e1 * t84) * t119) * t71) * t69) * t114 + (t84 * t116 + t62 * t120) * t56 + (-t61 + ((t86 - t91) * t82 * t78 * t71 * t69 + t137 * t112) * t62) * t84 * qJD(1)) * t58) * t93, 0, (t55 * t134 + t61 * t83) * t93 * t117 + ((t61 * t122 + (qJD(3) * t55 + t53) * t133) * t83 + (-t61 * t118 - (-t54 * t69 * t92 + t138 * t68 + (-qJD(3) * t68 - t121 * t69 + t131 * t57) * t66) * t114 + (t93 * t116 + t62 * t122) * t55 - ((-t54 + t121) * t68 + ((t66 * t92 - 0.1e1) * qJD(3) + (-t66 + t92) * t57) * t69) * t83 * t133) * t84) * t58, 0, 0; -0.2e1 * t108 * t83 * t140 + (qJD(3) * t102 + 0.2e1 * t101 * t83) * t74, 0, 0.2e1 * t84 * t125 * t140 + (t83 * t87 * t118 + qJD(1) * t102) * t74, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end