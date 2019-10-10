% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(10));
	t20 = sin(pkin(6));
	t19 = sin(pkin(10));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (27->18), div. (0->0), fcn. (42->8), ass. (0->14)
	t55 = sin(pkin(6));
	t60 = cos(qJ(2));
	t63 = t55 * t60;
	t58 = cos(pkin(6));
	t59 = sin(qJ(2));
	t62 = t58 * t59;
	t61 = t58 * t60;
	t57 = cos(pkin(10));
	t56 = cos(pkin(11));
	t54 = sin(pkin(10));
	t53 = sin(pkin(11));
	t52 = -t54 * t61 - t57 * t59;
	t51 = -t54 * t59 + t57 * t61;
	t1 = [0, t52 * t56, 0, 0, 0, 0; 0, t51 * t56, 0, 0, 0, 0; 0, t56 * t63, 0, 0, 0, 0; 0, -t52 * t53, 0, 0, 0, 0; 0, -t51 * t53, 0, 0, 0, 0; 0, -t53 * t63, 0, 0, 0, 0; 0, -t54 * t62 + t57 * t60, 0, 0, 0, 0; 0, t54 * t60 + t57 * t62, 0, 0, 0, 0; 0, t55 * t59, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (37->14), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->20)
	t72 = sin(pkin(10));
	t73 = sin(pkin(6));
	t83 = t72 * t73;
	t74 = cos(pkin(10));
	t82 = t73 * t74;
	t76 = sin(qJ(2));
	t81 = t73 * t76;
	t77 = cos(qJ(2));
	t80 = t73 * t77;
	t75 = cos(pkin(6));
	t79 = t75 * t76;
	t78 = t75 * t77;
	t71 = pkin(11) + qJ(4);
	t70 = cos(t71);
	t69 = sin(t71);
	t68 = -t72 * t79 + t74 * t77;
	t67 = -t72 * t78 - t74 * t76;
	t66 = t72 * t77 + t74 * t79;
	t65 = -t72 * t76 + t74 * t78;
	t1 = [0, t67 * t70, 0, -t68 * t69 + t70 * t83, 0, 0; 0, t65 * t70, 0, -t66 * t69 - t70 * t82, 0, 0; 0, t70 * t80, 0, -t69 * t81 + t75 * t70, 0, 0; 0, -t67 * t69, 0, -t68 * t70 - t69 * t83, 0, 0; 0, -t65 * t69, 0, -t66 * t70 + t69 * t82, 0, 0; 0, -t69 * t80, 0, -t75 * t69 - t70 * t81, 0, 0; 0, t68, 0, 0, 0, 0; 0, t66, 0, 0, 0, 0; 0, t81, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:21
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (37->14), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->20)
	t91 = sin(pkin(10));
	t92 = sin(pkin(6));
	t102 = t91 * t92;
	t93 = cos(pkin(10));
	t101 = t92 * t93;
	t95 = sin(qJ(2));
	t100 = t92 * t95;
	t96 = cos(qJ(2));
	t99 = t92 * t96;
	t94 = cos(pkin(6));
	t98 = t94 * t95;
	t97 = t94 * t96;
	t90 = pkin(11) + qJ(4);
	t89 = cos(t90);
	t88 = sin(t90);
	t87 = -t91 * t98 + t93 * t96;
	t86 = -t91 * t97 - t93 * t95;
	t85 = t91 * t96 + t93 * t98;
	t84 = -t91 * t95 + t93 * t97;
	t1 = [0, t87, 0, 0, 0, 0; 0, t85, 0, 0, 0, 0; 0, t100, 0, 0, 0, 0; 0, -t86 * t89, 0, -t89 * t102 + t87 * t88, 0, 0; 0, -t84 * t89, 0, t89 * t101 + t85 * t88, 0, 0; 0, -t89 * t99, 0, t88 * t100 - t94 * t89, 0, 0; 0, t86 * t88, 0, t88 * t102 + t87 * t89, 0, 0; 0, t84 * t88, 0, -t88 * t101 + t85 * t89, 0, 0; 0, t88 * t99, 0, t89 * t100 + t94 * t88, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:37:21
	% EndTime: 2019-10-09 21:37:22
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (90->28), mult. (163->65), div. (0->0), fcn. (238->10), ass. (0->31)
	t117 = pkin(11) + qJ(4);
	t115 = sin(t117);
	t122 = sin(qJ(6));
	t134 = t115 * t122;
	t124 = cos(qJ(6));
	t133 = t115 * t124;
	t118 = sin(pkin(10));
	t119 = sin(pkin(6));
	t132 = t118 * t119;
	t120 = cos(pkin(10));
	t131 = t119 * t120;
	t123 = sin(qJ(2));
	t130 = t119 * t123;
	t121 = cos(pkin(6));
	t129 = t121 * t123;
	t125 = cos(qJ(2));
	t128 = t121 * t125;
	t127 = t122 * t125;
	t126 = t124 * t125;
	t116 = cos(t117);
	t113 = -t118 * t129 + t120 * t125;
	t112 = t118 * t128 + t120 * t123;
	t111 = t118 * t125 + t120 * t129;
	t110 = t118 * t123 - t120 * t128;
	t109 = t121 * t115 + t116 * t130;
	t108 = t115 * t130 - t121 * t116;
	t107 = t113 * t116 + t115 * t132;
	t106 = t113 * t115 - t116 * t132;
	t105 = t111 * t116 - t115 * t131;
	t104 = t111 * t115 + t116 * t131;
	t1 = [0, -t112 * t134 + t113 * t124, 0, t107 * t122, 0, t106 * t124 - t112 * t122; 0, -t110 * t134 + t111 * t124, 0, t105 * t122, 0, t104 * t124 - t110 * t122; 0, (t115 * t127 + t123 * t124) * t119, 0, t109 * t122, 0, t108 * t124 + t119 * t127; 0, -t112 * t133 - t113 * t122, 0, t107 * t124, 0, -t106 * t122 - t112 * t124; 0, -t110 * t133 - t111 * t122, 0, t105 * t124, 0, -t104 * t122 - t110 * t124; 0, (t115 * t126 - t122 * t123) * t119, 0, t109 * t124, 0, -t108 * t122 + t119 * t126; 0, -t112 * t116, 0, -t106, 0, 0; 0, -t110 * t116, 0, -t104, 0, 0; 0, t119 * t125 * t116, 0, -t108, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end