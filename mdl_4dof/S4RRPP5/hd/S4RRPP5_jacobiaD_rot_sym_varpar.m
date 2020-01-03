% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RRPP5
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
%   Wie in S4RRPP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RRPP5_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPP5_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_jacobiaD_rot_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:00:46
	% EndTime: 2019-12-31 17:00:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:00:46
	% EndTime: 2019-12-31 17:00:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:00:46
	% EndTime: 2019-12-31 17:00:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:00:47
	% EndTime: 2019-12-31 17:00:47
	% DurationCPUTime: 0.51s
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
	t1 = [t87 * t79 * t102 + (qJD(2) * t97 - t117 * t120) * t66, t49, 0, 0; (t56 * t103 + (-t56 * t114 + (qJD(1) * t51 + t48) * t130) * t53) * t85 + (t57 * t103 * t51 + (-((t52 * t100 + t134 * t114 + t102) * t63 + (t101 * t124 - t131 + (t131 + (-t72 * t80 + t135) * t115) * t66) * t64) * t109 + (t84 * t111 - t57 * t114) * t51 + (-t56 + ((-t74 + t83) * t64 * t107 + t134 * t108) * t57) * t84 * qJD(1)) * t53) * t87, (t50 * t130 - t56 * t86) * t87 * t112 + ((-t56 * t117 + (-qJD(2) * t50 - t48) * t129) * t86 + (-t56 * t113 - (-t49 * t64 * t85 + t63 * t115 + t127 * t132 - t132 + (-qJD(2) * t63 - t116 * t64) * t62) * t109 + (t87 * t111 + t57 * t117) * t50 - ((t49 - t116) * t63 + ((0.1e1 - t127) * qJD(2) + (t62 - t85) * t52) * t64) * t86 * t129) * t84) * t53, 0, 0; t104 * t86 * t110 + (qJD(2) * t96 + 0.2e1 * t86 * t136) * t68, t75 * t110 * t119 + (-t75 * t86 * t113 + qJD(1) * t96) * t68, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:00:47
	% EndTime: 2019-12-31 17:00:47
	% DurationCPUTime: 0.51s
	% Computational Cost: add. (585->68), mult. (1839->161), div. (436->14), fcn. (2165->7), ass. (0->73)
	t87 = cos(qJ(1));
	t115 = qJD(1) * t87;
	t85 = sin(qJ(1));
	t78 = 0.1e1 / t85 ^ 2;
	t83 = t87 ^ 2;
	t122 = t78 * t83;
	t84 = sin(qJ(2));
	t72 = t84 ^ 2;
	t71 = t72 * t122 + 0.1e1;
	t76 = t85 ^ 2;
	t77 = 0.1e1 / t85;
	t95 = (-0.1e1 - 0.1e1 / t76 * t83) * t77 * t115;
	t114 = qJD(2) * t84;
	t86 = cos(qJ(2));
	t98 = t83 * t86 * t114;
	t135 = (t72 * t95 + t78 * t98) / t71 ^ 2;
	t74 = 0.1e1 / t84 ^ 2;
	t81 = t86 ^ 2;
	t123 = t74 * t81;
	t103 = 0.1e1 + t123;
	t133 = t103 * t85;
	t113 = qJD(2) * t85;
	t104 = t86 * t115;
	t70 = t76 * t123 + 0.1e1;
	t66 = 0.1e1 / t70;
	t73 = 0.1e1 / t84;
	t52 = ((t84 * t113 - t104) * t73 + t113 * t123) * t66;
	t132 = -t52 + t113;
	t119 = t85 * t86;
	t65 = atan2(-t119, t84);
	t63 = sin(t65);
	t106 = t63 * t119;
	t64 = cos(t65);
	t59 = t64 * t84 - t106;
	t56 = 0.1e1 / t59;
	t57 = 0.1e1 / t59 ^ 2;
	t131 = t66 - 0.1e1;
	t124 = t64 * t86;
	t48 = (-t52 * t85 + qJD(2)) * t124 + (t132 * t84 - t104) * t63;
	t130 = t48 * t56 * t57;
	t129 = t52 * t86;
	t128 = t57 * t86;
	t127 = t57 * t87;
	t80 = t86 * t81;
	t94 = qJD(2) * (-t86 - 0.1e1 / t72 * t80) * t73;
	t120 = t81 * t85;
	t99 = t115 * t120;
	t126 = (t74 * t99 + t76 * t94) / t70 ^ 2;
	t125 = t63 * t85;
	t121 = t81 * t83;
	t118 = t86 * t87;
	t117 = qJD(1) * t85;
	t116 = qJD(1) * t86;
	t112 = qJD(2) * t87;
	t55 = t57 * t121 + 0.1e1;
	t111 = 0.2e1 * (-t121 * t130 + (-t98 - t99) * t57) / t55 ^ 2;
	t110 = 0.2e1 * t130;
	t109 = -0.2e1 * t126;
	t108 = t57 * t118;
	t107 = t86 * t126;
	t105 = t73 * t120;
	t102 = 0.1e1 + t122;
	t101 = t86 * t111;
	t100 = t66 * t105;
	t97 = t103 * t87;
	t96 = t102 * t86;
	t68 = 0.1e1 / t71;
	t62 = t66 * t133;
	t53 = 0.1e1 / t55;
	t51 = (t131 * t86 * t63 + t64 * t100) * t87;
	t50 = t84 * t125 + t124 + (-t64 * t119 - t63 * t84) * t62;
	t49 = t109 * t133 + (qJD(1) * t97 + 0.2e1 * t85 * t94) * t66;
	t1 = [0.2e1 * t87 * t73 * t107 + (t73 * t85 * t116 + qJD(2) * t97) * t66, t49, 0, 0; (t56 * t101 + (t56 * t114 + (qJD(1) * t51 + t48) * t128) * t53) * t85 + (t57 * t101 * t51 + (-((-t52 * t100 - t131 * t114 - 0.2e1 * t107) * t63 + (t105 * t109 - t129 + (t129 + (-t74 * t80 - 0.2e1 * t86) * t113) * t66) * t64) * t108 + (t86 * t110 + t57 * t114) * t51 + (-t56 + ((t76 - t83) * t81 * t73 * t66 * t64 + t131 * t106) * t57) * t116) * t53) * t87, (t50 * t128 + t56 * t84) * t87 * t111 + ((t56 * t117 + (qJD(2) * t50 + t48) * t127) * t84 + (-t56 * t112 - (-t49 * t64 * t85 + t132 * t63 + (-qJD(2) * t63 - t115 * t64 + t125 * t52) * t62) * t108 + (t87 * t110 + t57 * t117) * t50 - ((-t49 + t115) * t63 + ((t62 * t85 - 0.1e1) * qJD(2) + (-t62 + t85) * t52) * t64) * t84 * t127) * t86) * t53, 0, 0; -0.2e1 * t102 * t84 * t135 + (qJD(2) * t96 + 0.2e1 * t84 * t95) * t68, 0.2e1 * t77 * t118 * t135 + (t77 * t84 * t112 + qJD(1) * t96) * t68, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end