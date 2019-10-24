% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRPRR6
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:27
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRPRR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:35
	% EndTime: 2019-10-24 10:27:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:35
	% EndTime: 2019-10-24 10:27:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:35
	% EndTime: 2019-10-24 10:27:35
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
	% StartTime: 2019-10-24 10:27:35
	% EndTime: 2019-10-24 10:27:35
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->7), mult. (27->18), div. (0->0), fcn. (42->8), ass. (0->14)
	t55 = sin(pkin(5));
	t60 = cos(qJ(2));
	t63 = t55 * t60;
	t58 = cos(pkin(5));
	t59 = sin(qJ(2));
	t62 = t58 * t59;
	t61 = t58 * t60;
	t57 = cos(pkin(9));
	t56 = cos(pkin(10));
	t54 = sin(pkin(9));
	t53 = sin(pkin(10));
	t52 = -t54 * t61 - t57 * t59;
	t51 = -t54 * t59 + t57 * t61;
	t1 = [0, t52 * t56, 0, 0, 0; 0, t51 * t56, 0, 0, 0; 0, t56 * t63, 0, 0, 0; 0, -t52 * t53, 0, 0, 0; 0, -t51 * t53, 0, 0, 0; 0, -t53 * t63, 0, 0, 0; 0, -t54 * t62 + t57 * t60, 0, 0, 0; 0, t54 * t60 + t57 * t62, 0, 0, 0; 0, t55 * t59, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:35
	% EndTime: 2019-10-24 10:27:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (37->14), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->20)
	t72 = sin(pkin(9));
	t73 = sin(pkin(5));
	t83 = t72 * t73;
	t74 = cos(pkin(9));
	t82 = t73 * t74;
	t76 = sin(qJ(2));
	t81 = t73 * t76;
	t77 = cos(qJ(2));
	t80 = t73 * t77;
	t75 = cos(pkin(5));
	t79 = t75 * t76;
	t78 = t75 * t77;
	t71 = pkin(10) + qJ(4);
	t70 = cos(t71);
	t69 = sin(t71);
	t68 = -t72 * t79 + t74 * t77;
	t67 = -t72 * t78 - t74 * t76;
	t66 = t72 * t77 + t74 * t79;
	t65 = -t72 * t76 + t74 * t78;
	t1 = [0, t67 * t70, 0, -t68 * t69 + t70 * t83, 0; 0, t65 * t70, 0, -t66 * t69 - t70 * t82, 0; 0, t70 * t80, 0, -t69 * t81 + t75 * t70, 0; 0, -t67 * t69, 0, -t68 * t70 - t69 * t83, 0; 0, -t65 * t69, 0, -t66 * t70 + t69 * t82, 0; 0, -t69 * t80, 0, -t75 * t69 - t70 * t81, 0; 0, t68, 0, 0, 0; 0, t66, 0, 0, 0; 0, t81, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:35
	% EndTime: 2019-10-24 10:27:35
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (93->28), mult. (163->65), div. (0->0), fcn. (238->10), ass. (0->31)
	t115 = pkin(10) + qJ(4);
	t114 = cos(t115);
	t120 = sin(qJ(5));
	t132 = t114 * t120;
	t122 = cos(qJ(5));
	t131 = t114 * t122;
	t116 = sin(pkin(9));
	t117 = sin(pkin(5));
	t130 = t116 * t117;
	t118 = cos(pkin(9));
	t129 = t117 * t118;
	t121 = sin(qJ(2));
	t128 = t117 * t121;
	t119 = cos(pkin(5));
	t127 = t119 * t121;
	t123 = cos(qJ(2));
	t126 = t119 * t123;
	t125 = t120 * t123;
	t124 = t122 * t123;
	t113 = sin(t115);
	t111 = -t116 * t127 + t118 * t123;
	t110 = t116 * t126 + t118 * t121;
	t109 = t116 * t123 + t118 * t127;
	t108 = t116 * t121 - t118 * t126;
	t107 = t119 * t113 + t114 * t128;
	t106 = -t113 * t128 + t119 * t114;
	t105 = t111 * t114 + t113 * t130;
	t104 = -t111 * t113 + t114 * t130;
	t103 = t109 * t114 - t113 * t129;
	t102 = -t109 * t113 - t114 * t129;
	t1 = [0, -t110 * t131 + t111 * t120, 0, t104 * t122, -t105 * t120 + t110 * t122; 0, -t108 * t131 + t109 * t120, 0, t102 * t122, -t103 * t120 + t108 * t122; 0, (t114 * t124 + t120 * t121) * t117, 0, t106 * t122, -t107 * t120 - t117 * t124; 0, t110 * t132 + t111 * t122, 0, -t104 * t120, -t105 * t122 - t110 * t120; 0, t108 * t132 + t109 * t122, 0, -t102 * t120, -t103 * t122 - t108 * t120; 0, (-t114 * t125 + t121 * t122) * t117, 0, -t106 * t120, -t107 * t122 + t117 * t125; 0, -t110 * t113, 0, t105, 0; 0, -t108 * t113, 0, t103, 0; 0, t117 * t123 * t113, 0, t107, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,5);
end