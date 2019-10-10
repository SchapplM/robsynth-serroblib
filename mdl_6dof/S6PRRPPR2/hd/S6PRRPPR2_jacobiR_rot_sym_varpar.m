% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:09
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR2_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR2_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
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
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (19->13), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->19)
	t66 = sin(pkin(6));
	t69 = sin(qJ(3));
	t78 = t66 * t69;
	t70 = sin(qJ(2));
	t77 = t66 * t70;
	t71 = cos(qJ(3));
	t76 = t66 * t71;
	t72 = cos(qJ(2));
	t75 = t66 * t72;
	t68 = cos(pkin(6));
	t74 = t68 * t70;
	t73 = t68 * t72;
	t67 = cos(pkin(10));
	t65 = sin(pkin(10));
	t64 = -t65 * t74 + t67 * t72;
	t63 = -t65 * t73 - t67 * t70;
	t62 = t65 * t72 + t67 * t74;
	t61 = -t65 * t70 + t67 * t73;
	t1 = [0, t63 * t71, -t64 * t69 + t65 * t76, 0, 0, 0; 0, t61 * t71, -t62 * t69 - t67 * t76, 0, 0, 0; 0, t71 * t75, t68 * t71 - t69 * t77, 0, 0, 0; 0, -t63 * t69, -t64 * t71 - t65 * t78, 0, 0, 0; 0, -t61 * t69, -t62 * t71 + t67 * t78, 0, 0, 0; 0, -t69 * t75, -t68 * t69 - t70 * t76, 0, 0, 0; 0, t64, 0, 0, 0, 0; 0, t62, 0, 0, 0, 0; 0, t77, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (37->14), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->20)
	t75 = sin(pkin(10));
	t76 = sin(pkin(6));
	t86 = t75 * t76;
	t77 = cos(pkin(10));
	t85 = t76 * t77;
	t79 = sin(qJ(2));
	t84 = t76 * t79;
	t80 = cos(qJ(2));
	t83 = t76 * t80;
	t78 = cos(pkin(6));
	t82 = t78 * t79;
	t81 = t78 * t80;
	t74 = qJ(3) + pkin(11);
	t73 = cos(t74);
	t72 = sin(t74);
	t71 = -t75 * t82 + t77 * t80;
	t70 = -t75 * t81 - t77 * t79;
	t69 = t75 * t80 + t77 * t82;
	t68 = -t75 * t79 + t77 * t81;
	t1 = [0, t70 * t73, -t71 * t72 + t73 * t86, 0, 0, 0; 0, t68 * t73, -t69 * t72 - t73 * t85, 0, 0, 0; 0, t73 * t83, -t72 * t84 + t73 * t78, 0, 0, 0; 0, -t70 * t72, -t71 * t73 - t72 * t86, 0, 0, 0; 0, -t68 * t72, -t69 * t73 + t72 * t85, 0, 0, 0; 0, -t72 * t83, -t72 * t78 - t73 * t84, 0, 0, 0; 0, t71, 0, 0, 0, 0; 0, t69, 0, 0, 0, 0; 0, t84, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (37->14), mult. (57->32), div. (0->0), fcn. (88->8), ass. (0->20)
	t94 = sin(pkin(10));
	t95 = sin(pkin(6));
	t105 = t94 * t95;
	t96 = cos(pkin(10));
	t104 = t95 * t96;
	t98 = sin(qJ(2));
	t103 = t95 * t98;
	t99 = cos(qJ(2));
	t102 = t95 * t99;
	t97 = cos(pkin(6));
	t101 = t97 * t98;
	t100 = t97 * t99;
	t93 = qJ(3) + pkin(11);
	t92 = cos(t93);
	t91 = sin(t93);
	t90 = -t94 * t101 + t96 * t99;
	t89 = -t94 * t100 - t96 * t98;
	t88 = t96 * t101 + t94 * t99;
	t87 = t96 * t100 - t94 * t98;
	t1 = [0, t90, 0, 0, 0, 0; 0, t88, 0, 0, 0, 0; 0, t103, 0, 0, 0, 0; 0, -t89 * t92, -t92 * t105 + t90 * t91, 0, 0, 0; 0, -t87 * t92, t92 * t104 + t88 * t91, 0, 0, 0; 0, -t92 * t102, t91 * t103 - t97 * t92, 0, 0, 0; 0, t89 * t91, t91 * t105 + t90 * t92, 0, 0, 0; 0, t87 * t91, -t91 * t104 + t88 * t92, 0, 0, 0; 0, t91 * t102, t92 * t103 + t97 * t91, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (90->28), mult. (163->65), div. (0->0), fcn. (238->10), ass. (0->31)
	t121 = qJ(3) + pkin(11);
	t119 = sin(t121);
	t126 = sin(qJ(6));
	t138 = t119 * t126;
	t128 = cos(qJ(6));
	t137 = t119 * t128;
	t122 = sin(pkin(10));
	t123 = sin(pkin(6));
	t136 = t122 * t123;
	t124 = cos(pkin(10));
	t135 = t123 * t124;
	t127 = sin(qJ(2));
	t134 = t123 * t127;
	t125 = cos(pkin(6));
	t133 = t125 * t127;
	t129 = cos(qJ(2));
	t132 = t125 * t129;
	t131 = t126 * t129;
	t130 = t128 * t129;
	t120 = cos(t121);
	t117 = -t122 * t133 + t124 * t129;
	t116 = t122 * t132 + t124 * t127;
	t115 = t122 * t129 + t124 * t133;
	t114 = t122 * t127 - t124 * t132;
	t113 = t125 * t119 + t120 * t134;
	t112 = t119 * t134 - t125 * t120;
	t111 = t117 * t120 + t119 * t136;
	t110 = t117 * t119 - t120 * t136;
	t109 = t115 * t120 - t119 * t135;
	t108 = t115 * t119 + t120 * t135;
	t1 = [0, -t116 * t138 + t117 * t128, t111 * t126, 0, 0, t110 * t128 - t116 * t126; 0, -t114 * t138 + t115 * t128, t109 * t126, 0, 0, t108 * t128 - t114 * t126; 0, (t119 * t131 + t127 * t128) * t123, t113 * t126, 0, 0, t112 * t128 + t123 * t131; 0, -t116 * t137 - t117 * t126, t111 * t128, 0, 0, -t110 * t126 - t116 * t128; 0, -t114 * t137 - t115 * t126, t109 * t128, 0, 0, -t108 * t126 - t114 * t128; 0, (t119 * t130 - t126 * t127) * t123, t113 * t128, 0, 0, -t112 * t126 + t123 * t130; 0, -t116 * t120, -t110, 0, 0, 0; 0, -t114 * t120, -t108, 0, 0, 0; 0, t123 * t129 * t120, -t112, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end