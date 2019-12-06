% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRPR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:33:48
	% EndTime: 2019-12-05 16:33:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:33:48
	% EndTime: 2019-12-05 16:33:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:33:48
	% EndTime: 2019-12-05 16:33:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(5));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(9));
	t20 = sin(pkin(5));
	t19 = sin(pkin(9));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0; 0, t20 * t24, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0; 0, -t20 * t23, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:33:48
	% EndTime: 2019-12-05 16:33:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (19->13), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t66 = sin(pkin(5));
	t69 = sin(qJ(3));
	t78 = t66 * t69;
	t70 = sin(qJ(2));
	t77 = t66 * t70;
	t71 = cos(qJ(3));
	t76 = t66 * t71;
	t72 = cos(qJ(2));
	t75 = t66 * t72;
	t68 = cos(pkin(5));
	t74 = t68 * t70;
	t73 = t68 * t72;
	t67 = cos(pkin(9));
	t65 = sin(pkin(9));
	t64 = -t65 * t74 + t67 * t72;
	t63 = -t65 * t73 - t67 * t70;
	t62 = t65 * t72 + t67 * t74;
	t61 = -t65 * t70 + t67 * t73;
	t1 = [0, t63 * t71, -t64 * t69 + t65 * t76, 0, 0; 0, t61 * t71, -t62 * t69 - t67 * t76, 0, 0; 0, t71 * t75, t68 * t71 - t69 * t77, 0, 0; 0, -t63 * t69, -t64 * t71 - t65 * t78, 0, 0; 0, -t61 * t69, -t62 * t71 + t67 * t78, 0, 0; 0, -t69 * t75, -t68 * t69 - t70 * t76, 0, 0; 0, t64, 0, 0, 0; 0, t62, 0, 0, 0; 0, t77, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:33:49
	% EndTime: 2019-12-05 16:33:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (34->19), mult. (107->48), div. (0->0), fcn. (156->10), ass. (0->27)
	t100 = sin(qJ(3));
	t96 = sin(pkin(5));
	t112 = t100 * t96;
	t102 = cos(qJ(3));
	t94 = sin(pkin(10));
	t111 = t102 * t94;
	t110 = t102 * t96;
	t97 = cos(pkin(10));
	t109 = t102 * t97;
	t101 = sin(qJ(2));
	t95 = sin(pkin(9));
	t108 = t95 * t101;
	t103 = cos(qJ(2));
	t107 = t95 * t103;
	t98 = cos(pkin(9));
	t106 = t98 * t101;
	t105 = t98 * t103;
	t104 = t102 * t103;
	t99 = cos(pkin(5));
	t93 = -t101 * t112 + t99 * t102;
	t92 = -t99 * t108 + t105;
	t91 = -t99 * t107 - t106;
	t90 = t99 * t106 + t107;
	t89 = t99 * t105 - t108;
	t88 = -t92 * t100 + t95 * t110;
	t87 = -t90 * t100 - t98 * t110;
	t1 = [0, t91 * t109 + t92 * t94, t88 * t97, 0, 0; 0, t89 * t109 + t90 * t94, t87 * t97, 0, 0; 0, (t101 * t94 + t97 * t104) * t96, t93 * t97, 0, 0; 0, -t91 * t111 + t92 * t97, -t88 * t94, 0, 0; 0, -t89 * t111 + t90 * t97, -t87 * t94, 0, 0; 0, (t101 * t97 - t94 * t104) * t96, -t93 * t94, 0, 0; 0, t91 * t100, t92 * t102 + t95 * t112, 0, 0; 0, t89 * t100, t90 * t102 - t98 * t112, 0, 0; 0, t103 * t112, t99 * t100 + t101 * t110, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:33:49
	% EndTime: 2019-12-05 16:33:49
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (84->28), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->30)
	t115 = pkin(10) + qJ(5);
	t113 = sin(t115);
	t122 = cos(qJ(3));
	t131 = t113 * t122;
	t114 = cos(t115);
	t130 = t114 * t122;
	t117 = sin(pkin(5));
	t120 = sin(qJ(3));
	t129 = t117 * t120;
	t128 = t117 * t122;
	t123 = cos(qJ(2));
	t127 = t117 * t123;
	t119 = cos(pkin(5));
	t121 = sin(qJ(2));
	t126 = t119 * t121;
	t125 = t119 * t123;
	t124 = t122 * t123;
	t118 = cos(pkin(9));
	t116 = sin(pkin(9));
	t111 = t119 * t120 + t121 * t128;
	t110 = t119 * t122 - t121 * t129;
	t109 = -t116 * t126 + t118 * t123;
	t108 = t116 * t125 + t118 * t121;
	t107 = t116 * t123 + t118 * t126;
	t106 = t116 * t121 - t118 * t125;
	t105 = t109 * t122 + t116 * t129;
	t104 = -t109 * t120 + t116 * t128;
	t103 = t107 * t122 - t118 * t129;
	t102 = -t107 * t120 - t118 * t128;
	t1 = [0, -t108 * t130 + t109 * t113, t104 * t114, 0, -t105 * t113 + t108 * t114; 0, -t106 * t130 + t107 * t113, t102 * t114, 0, -t103 * t113 + t106 * t114; 0, (t113 * t121 + t114 * t124) * t117, t110 * t114, 0, -t111 * t113 - t114 * t127; 0, t108 * t131 + t109 * t114, -t104 * t113, 0, -t105 * t114 - t108 * t113; 0, t106 * t131 + t107 * t114, -t102 * t113, 0, -t103 * t114 - t106 * t113; 0, (-t113 * t124 + t114 * t121) * t117, -t110 * t113, 0, -t111 * t114 + t113 * t127; 0, -t108 * t120, t105, 0, 0; 0, -t106 * t120, t103, 0, 0; 0, t120 * t127, t111, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end