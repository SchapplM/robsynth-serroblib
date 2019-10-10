% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4RPPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RPPP1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPP1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:02
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (137->30), mult. (614->95), div. (108->12), fcn. (792->9), ass. (0->49)
	t86 = sin(pkin(4));
	t79 = t86 ^ 2;
	t88 = cos(pkin(4));
	t81 = 0.1e1 / t88 ^ 2;
	t90 = cos(qJ(1));
	t84 = t90 ^ 2;
	t77 = t79 * t81 * t84 + 0.1e1;
	t89 = sin(qJ(1));
	t83 = t89 ^ 2;
	t108 = 0.1e1 / t77 ^ 2 * t83;
	t112 = t108 * t81;
	t103 = t90 * t86;
	t76 = atan2(t103, t88);
	t72 = sin(t76);
	t73 = cos(t76);
	t58 = t72 * t103 + t73 * t88;
	t55 = 0.1e1 / t58;
	t85 = sin(pkin(6));
	t105 = t89 * t85;
	t87 = cos(pkin(6));
	t99 = t88 * t105 - t87 * t90;
	t65 = 0.1e1 / t99;
	t80 = 0.1e1 / t88;
	t56 = 0.1e1 / t58 ^ 2;
	t66 = 0.1e1 / t99 ^ 2;
	t111 = t56 * t89;
	t104 = t89 * t87;
	t70 = t88 * t104 + t85 * t90;
	t110 = t66 * t70;
	t106 = t88 * t90;
	t69 = -t85 * t106 - t104;
	t109 = t69 * t70;
	t107 = t79 * t80;
	t102 = qJD(1) * t90;
	t74 = 0.1e1 / t77;
	t101 = (t74 - 0.1e1) * t86;
	t100 = -0.2e1 * t80 * t112;
	t68 = t87 * t106 - t105;
	t51 = (-t73 * t74 * t90 * t107 + t72 * t101) * t89;
	t78 = t86 * t79;
	t67 = t65 * t66;
	t64 = t70 ^ 2;
	t63 = t69 * qJD(1);
	t62 = t68 * qJD(1);
	t61 = t64 * t66 + 0.1e1;
	t57 = t55 * t56;
	t54 = t56 * t79 * t83 + 0.1e1;
	t50 = qJD(1) * t51;
	t1 = [(-t74 * t80 * t86 + t78 * t100) * t102, 0, 0, 0; (0.2e1 * (t51 * t111 - t55 * t90) / t54 ^ 2 * (-t50 * t57 * t83 + t102 * t111) * t79 + ((0.2e1 * t51 * t57 * t89 - t56 * t90) * t50 + (-t89 * t55 + ((-t51 + (-t78 * t112 - t101) * t89 * t72) * t90 - (t84 * t79 ^ 2 * t100 + (-t108 + (0.2e1 * t83 - t84) * t74) * t107) * t89 * t73) * t56) * qJD(1)) / t54) * t86, 0, 0, 0; 0.2e1 * (t66 * t109 + t65 * t68) / t61 ^ 2 * (t63 * t64 * t67 + t62 * t110) + (-t69 * t62 * t66 + (-0.2e1 * t67 * t109 - t66 * t68) * t63 + (-t99 * t110 + t70 * t65) * qJD(1)) / t61, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:02
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (361->28), mult. (1032->88), div. (172->14), fcn. (1422->9), ass. (0->50)
	t100 = cos(qJ(1));
	t120 = cos(pkin(4));
	t99 = sin(qJ(1));
	t116 = t99 * t120;
	t96 = sin(pkin(6));
	t98 = cos(pkin(6));
	t110 = -t100 * t98 + t96 * t116;
	t79 = t110 ^ 2;
	t97 = sin(pkin(4));
	t88 = 0.1e1 / t97 ^ 2;
	t94 = 0.1e1 / t99 ^ 2;
	t74 = t79 * t94 * t88 + 0.1e1;
	t135 = -0.2e1 / t74;
	t93 = 0.1e1 / t99;
	t133 = qJD(1) * t93 * t94;
	t115 = t100 * t120;
	t82 = -t96 * t115 - t99 * t98;
	t76 = t82 * qJD(1);
	t134 = -0.2e1 * (-t100 * t79 * t133 - t110 * t76 * t94) * t88 / t74 ^ 2;
	t87 = 0.1e1 / t97;
	t125 = t87 / t98;
	t124 = t88 / t98 ^ 2;
	t80 = -t98 * t115 + t99 * t96;
	t122 = t97 * t98;
	t68 = atan2(-t80, -t122);
	t66 = sin(t68);
	t67 = cos(t68);
	t130 = t80 ^ 2;
	t73 = t130 * t124 + 0.1e1;
	t69 = 0.1e1 / t73;
	t131 = (t67 * t80 * t125 - t66) * t69 + t66;
	t64 = -t67 * t122 - t66 * t80;
	t61 = 0.1e1 / t64;
	t62 = 0.1e1 / t64 ^ 2;
	t83 = t100 * t96 + t98 * t116;
	t77 = t83 * qJD(1);
	t56 = t131 * t77;
	t129 = t56 * t61 * t62;
	t70 = 0.1e1 / t73 ^ 2;
	t128 = t70 * t80;
	t75 = t80 * qJD(1);
	t127 = t75 * t62;
	t126 = t77 * t83;
	t123 = t124 * t125;
	t119 = t69 * t125;
	t118 = t100 * t87 * t94;
	t78 = t83 ^ 2;
	t60 = t78 * t62 + 0.1e1;
	t57 = t131 * t83;
	t1 = [-0.2e1 * t123 * t126 * t128 - t75 * t119, 0, 0, 0; 0.2e1 * (-t57 * t62 * t83 + t61 * t80) / t60 ^ 2 * (-t83 * t127 + t78 * t129) + (-t77 * t61 + (-t80 * t56 - t57 * t75) * t62 + (0.2e1 * t57 * t129 - t131 * t127 - (-t66 * t124 * t128 + (-0.2e1 * t119 + (0.2e1 * t130 * t123 + t125) * t70) * t67) * t62 * t126) * t83) / t60, 0, 0, 0; t82 * t93 * t87 * t134 + t76 * t118 * t135 + (t100 ^ 2 * t87 * t133 * t135 + t118 * t134) * t110, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:02
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (295->28), mult. (1032->91), div. (172->14), fcn. (1422->9), ass. (0->52)
	t100 = cos(qJ(1));
	t121 = cos(pkin(4));
	t99 = sin(qJ(1));
	t116 = t99 * t121;
	t96 = sin(pkin(6));
	t98 = cos(pkin(6));
	t110 = t100 * t96 + t98 * t116;
	t78 = t110 ^ 2;
	t97 = sin(pkin(4));
	t91 = 0.1e1 / t97 ^ 2;
	t94 = 0.1e1 / t99 ^ 2;
	t74 = t78 * t94 * t91 + 0.1e1;
	t133 = -0.2e1 / t74;
	t120 = qJD(1) * t100;
	t115 = t100 * t121;
	t80 = -t98 * t115 + t99 * t96;
	t75 = t80 * qJD(1);
	t93 = 0.1e1 / t99;
	t95 = t93 * t94;
	t132 = -0.2e1 * (-t110 * t75 * t94 - t78 * t95 * t120) * t91 / t74 ^ 2;
	t90 = 0.1e1 / t97;
	t125 = 0.1e1 / t96 * t90;
	t124 = 0.1e1 / t96 ^ 2 * t91;
	t122 = t97 * t96;
	t81 = t96 * t115 + t99 * t98;
	t70 = atan2(-t81, t122);
	t66 = sin(t70);
	t67 = cos(t70);
	t130 = t81 ^ 2;
	t73 = t130 * t124 + 0.1e1;
	t68 = 0.1e1 / t73;
	t109 = (t67 * t81 * t125 + t66) * t68 - t66;
	t64 = t67 * t122 - t66 * t81;
	t61 = 0.1e1 / t64;
	t62 = 0.1e1 / t64 ^ 2;
	t112 = t96 * t116;
	t77 = -qJD(1) * t112 + t98 * t120;
	t56 = t109 * t77;
	t129 = t56 * t61 * t62;
	t69 = 0.1e1 / t73 ^ 2;
	t128 = t69 * t81;
	t76 = t81 * qJD(1);
	t127 = t76 * t62;
	t84 = t100 * t98 - t112;
	t126 = t77 * t84;
	t123 = t124 * t125;
	t119 = t68 * t125;
	t118 = t100 * t90 * t94;
	t79 = t84 ^ 2;
	t60 = t79 * t62 + 0.1e1;
	t57 = t109 * t84;
	t1 = [0.2e1 * t123 * t126 * t128 + t76 * t119, 0, 0, 0; 0.2e1 * (t57 * t62 * t84 + t61 * t81) / t60 ^ 2 * (-t84 * t127 - t79 * t129) + (-t77 * t61 + (t81 * t56 + t57 * t76) * t62 + (0.2e1 * t57 * t129 + t109 * t127 - (-t66 * t124 * t128 + (0.2e1 * t119 + (-0.2e1 * t130 * t123 - t125) * t69) * t67) * t62 * t126) * t84) / t60, 0, 0, 0; t80 * t93 * t90 * t132 + t75 * t118 * t133 + (t100 ^ 2 * t90 * t95 * qJD(1) * t133 + t118 * t132) * t110, 0, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end