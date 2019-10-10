% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiR_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t50 = sin(qJ(2));
	t51 = sin(qJ(1));
	t57 = t51 * t50;
	t52 = cos(qJ(2));
	t56 = t51 * t52;
	t53 = cos(qJ(1));
	t55 = t53 * t50;
	t54 = t53 * t52;
	t49 = cos(pkin(6));
	t48 = sin(pkin(6));
	t47 = -t49 * t57 + t54;
	t46 = -t49 * t56 - t55;
	t45 = -t49 * t55 - t56;
	t44 = -t49 * t54 + t57;
	t1 = [t45, t46, 0, 0, 0, 0; t47, -t44, 0, 0, 0, 0; 0, t48 * t52, 0, 0, 0, 0; t44, -t47, 0, 0, 0, 0; t46, t45, 0, 0, 0, 0; 0, -t48 * t50, 0, 0, 0, 0; t53 * t48, 0, 0, 0, 0, 0; t51 * t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (8->6), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t61 = sin(qJ(2));
	t62 = sin(qJ(1));
	t68 = t62 * t61;
	t63 = cos(qJ(2));
	t67 = t62 * t63;
	t64 = cos(qJ(1));
	t66 = t64 * t61;
	t65 = t64 * t63;
	t60 = cos(pkin(6));
	t59 = sin(pkin(6));
	t58 = -t60 * t68 + t65;
	t57 = t60 * t67 + t66;
	t56 = t60 * t66 + t67;
	t55 = t60 * t65 - t68;
	t1 = [-t56, -t57, 0, 0, 0, 0; t58, t55, 0, 0, 0, 0; 0, t59 * t63, 0, 0, 0, 0; t64 * t59, 0, 0, 0, 0, 0; t62 * t59, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t55, t58, 0, 0, 0, 0; t57, t56, 0, 0, 0, 0; 0, t59 * t61, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t56 = sin(qJ(2));
	t57 = sin(qJ(1));
	t63 = t57 * t56;
	t58 = cos(qJ(2));
	t62 = t57 * t58;
	t59 = cos(qJ(1));
	t61 = t59 * t56;
	t60 = t59 * t58;
	t55 = cos(pkin(6));
	t54 = sin(pkin(6));
	t53 = -t55 * t63 + t60;
	t52 = t55 * t62 + t61;
	t51 = t55 * t61 + t62;
	t50 = -t55 * t60 + t63;
	t1 = [-t50, t53, 0, 0, 0, 0; t52, t51, 0, 0, 0, 0; 0, t54 * t56, 0, 0, 0, 0; t51, t52, 0, 0, 0, 0; -t53, t50, 0, 0, 0, 0; 0, -t54 * t58, 0, 0, 0, 0; -t59 * t54, 0, 0, 0, 0, 0; -t57 * t54, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (29->18), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->25)
	t82 = sin(pkin(6));
	t84 = sin(qJ(5));
	t99 = t82 * t84;
	t87 = cos(qJ(5));
	t98 = t82 * t87;
	t88 = cos(qJ(2));
	t97 = t82 * t88;
	t89 = cos(qJ(1));
	t96 = t82 * t89;
	t85 = sin(qJ(2));
	t86 = sin(qJ(1));
	t95 = t86 * t85;
	t94 = t86 * t88;
	t93 = t89 * t85;
	t92 = t89 * t88;
	t83 = cos(pkin(6));
	t77 = -t83 * t92 + t95;
	t91 = -t77 * t84 + t87 * t96;
	t90 = -t77 * t87 - t84 * t96;
	t80 = -t83 * t95 + t92;
	t79 = t83 * t94 + t93;
	t78 = t83 * t93 + t94;
	t76 = t79 * t87 - t86 * t99;
	t75 = -t79 * t84 - t86 * t98;
	t1 = [t90, t80 * t87, 0, 0, t75, 0; t76, t78 * t87, 0, 0, t91, 0; 0, t85 * t98, 0, 0, -t83 * t87 + t84 * t97, 0; -t91, -t80 * t84, 0, 0, -t76, 0; t75, -t78 * t84, 0, 0, t90, 0; 0, -t85 * t99, 0, 0, t83 * t84 + t87 * t97, 0; -t78, -t79, 0, 0, 0, 0; t80, -t77, 0, 0, 0, 0; 0, t97, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:41
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (74->27), mult. (219->61), div. (0->0), fcn. (320->10), ass. (0->35)
	t127 = cos(pkin(6));
	t134 = cos(qJ(2));
	t135 = cos(qJ(1));
	t136 = t135 * t134;
	t130 = sin(qJ(2));
	t131 = sin(qJ(1));
	t140 = t131 * t130;
	t120 = -t127 * t136 + t140;
	t129 = sin(qJ(5));
	t133 = cos(qJ(5));
	t126 = sin(pkin(6));
	t142 = t126 * t135;
	t114 = t120 * t133 + t129 * t142;
	t137 = t135 * t130;
	t139 = t131 * t134;
	t121 = t127 * t137 + t139;
	t128 = sin(qJ(6));
	t132 = cos(qJ(6));
	t149 = t114 * t128 - t121 * t132;
	t148 = -t114 * t132 - t121 * t128;
	t145 = t126 * t130;
	t144 = t126 * t131;
	t143 = t126 * t134;
	t141 = t128 * t133;
	t138 = t132 * t133;
	t113 = -t120 * t129 + t133 * t142;
	t123 = -t127 * t140 + t136;
	t122 = t127 * t139 + t137;
	t119 = -t127 * t129 - t133 * t143;
	t118 = -t127 * t133 + t129 * t143;
	t117 = t122 * t133 - t129 * t144;
	t116 = t122 * t129 + t133 * t144;
	t112 = t117 * t132 + t123 * t128;
	t111 = -t117 * t128 + t123 * t132;
	t1 = [t148, -t122 * t128 + t123 * t138, 0, 0, -t116 * t132, t111; t112, -t120 * t128 + t121 * t138, 0, 0, t113 * t132, -t149; 0, (t128 * t134 + t130 * t138) * t126, 0, 0, t118 * t132, -t119 * t128 + t132 * t145; t149, -t122 * t132 - t123 * t141, 0, 0, t116 * t128, -t112; t111, -t120 * t132 - t121 * t141, 0, 0, -t113 * t128, t148; 0, (-t130 * t141 + t132 * t134) * t126, 0, 0, -t118 * t128, -t119 * t132 - t128 * t145; t113, t123 * t129, 0, 0, t117, 0; t116, t121 * t129, 0, 0, t114, 0; 0, t129 * t145, 0, 0, t119, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end