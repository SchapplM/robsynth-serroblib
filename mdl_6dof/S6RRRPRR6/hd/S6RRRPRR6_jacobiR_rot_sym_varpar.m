% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:02
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
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
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (8->4), div. (0->0), fcn. (20->4), ass. (0->9)
	t10 = sin(qJ(1));
	t9 = sin(qJ(2));
	t16 = t10 * t9;
	t12 = cos(qJ(1));
	t15 = t12 * t9;
	t11 = cos(qJ(2));
	t14 = t10 * t11;
	t13 = t12 * t11;
	t1 = [-t14, -t15, 0, 0, 0, 0; t13, -t16, 0, 0, 0, 0; 0, t11, 0, 0, 0, 0; t16, -t13, 0, 0, 0, 0; -t15, -t14, 0, 0, 0, 0; 0, -t9, 0, 0, 0, 0; t12, 0, 0, 0, 0, 0; t10, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (14->12), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->18)
	t60 = sin(qJ(2));
	t61 = sin(qJ(1));
	t71 = t61 * t60;
	t62 = cos(qJ(3));
	t70 = t61 * t62;
	t59 = sin(qJ(3));
	t63 = cos(qJ(2));
	t69 = t63 * t59;
	t68 = t63 * t62;
	t64 = cos(qJ(1));
	t67 = t64 * t60;
	t66 = t64 * t62;
	t65 = t64 * t63;
	t58 = t61 * t59 + t62 * t65;
	t57 = -t59 * t65 + t70;
	t56 = t64 * t59 - t61 * t68;
	t55 = t61 * t69 + t66;
	t1 = [t56, -t60 * t66, t57, 0, 0, 0; t58, -t60 * t70, -t55, 0, 0, 0; 0, t68, -t60 * t59, 0, 0, 0; t55, t59 * t67, -t58, 0, 0, 0; t57, t59 * t71, t56, 0, 0, 0; 0, -t69, -t60 * t62, 0, 0, 0; -t71, t65, 0, 0, 0, 0; t67, t61 * t63, 0, 0, 0, 0; 0, t60, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (38->13), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->17)
	t70 = sin(qJ(2));
	t71 = sin(qJ(1));
	t78 = t71 * t70;
	t69 = qJ(3) + pkin(11);
	t67 = sin(t69);
	t72 = cos(qJ(2));
	t77 = t72 * t67;
	t68 = cos(t69);
	t76 = t72 * t68;
	t73 = cos(qJ(1));
	t75 = t73 * t70;
	t74 = t73 * t72;
	t66 = t71 * t67 + t68 * t74;
	t65 = -t67 * t74 + t71 * t68;
	t64 = t73 * t67 - t71 * t76;
	t63 = t73 * t68 + t71 * t77;
	t1 = [t64, -t68 * t75, t65, 0, 0, 0; t66, -t68 * t78, -t63, 0, 0, 0; 0, t76, -t70 * t67, 0, 0, 0; t63, t67 * t75, -t66, 0, 0, 0; t65, t67 * t78, t64, 0, 0, 0; 0, -t77, -t70 * t68, 0, 0, 0; -t78, t74, 0, 0, 0, 0; t75, t71 * t72, 0, 0, 0, 0; 0, t70, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (88->17), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->19)
	t83 = qJ(3) + pkin(11) + qJ(5);
	t81 = sin(t83);
	t84 = sin(qJ(2));
	t94 = t84 * t81;
	t82 = cos(t83);
	t93 = t84 * t82;
	t85 = sin(qJ(1));
	t92 = t85 * t84;
	t86 = cos(qJ(2));
	t91 = t86 * t81;
	t90 = t86 * t82;
	t87 = cos(qJ(1));
	t89 = t87 * t84;
	t88 = t87 * t86;
	t80 = t85 * t81 + t82 * t88;
	t79 = -t81 * t88 + t85 * t82;
	t78 = t87 * t81 - t85 * t90;
	t77 = t87 * t82 + t85 * t91;
	t1 = [t78, -t82 * t89, t79, 0, t79, 0; t80, -t82 * t92, -t77, 0, -t77, 0; 0, t90, -t94, 0, -t94, 0; t77, t81 * t89, -t80, 0, -t80, 0; t79, t81 * t92, t78, 0, t78, 0; 0, -t91, -t93, 0, -t93, 0; -t92, t88, 0, 0, 0, 0; t89, t85 * t86, 0, 0, 0, 0; 0, t84, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:02:21
	% EndTime: 2019-10-10 12:02:21
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (158->21), mult. (68->20), div. (0->0), fcn. (117->6), ass. (0->19)
	t89 = qJ(3) + pkin(11) + qJ(5) + qJ(6);
	t87 = sin(t89);
	t90 = sin(qJ(2));
	t100 = t90 * t87;
	t88 = cos(t89);
	t99 = t90 * t88;
	t91 = sin(qJ(1));
	t98 = t91 * t90;
	t92 = cos(qJ(2));
	t97 = t92 * t87;
	t96 = t92 * t88;
	t93 = cos(qJ(1));
	t95 = t93 * t90;
	t94 = t93 * t92;
	t86 = t91 * t87 + t88 * t94;
	t85 = -t87 * t94 + t91 * t88;
	t84 = t93 * t87 - t91 * t96;
	t83 = t93 * t88 + t91 * t97;
	t1 = [t84, -t88 * t95, t85, 0, t85, t85; t86, -t88 * t98, -t83, 0, -t83, -t83; 0, t96, -t100, 0, -t100, -t100; t83, t87 * t95, -t86, 0, -t86, -t86; t85, t87 * t98, t84, 0, t84, t84; 0, -t97, -t99, 0, -t99, -t99; -t98, t94, 0, 0, 0, 0; t95, t91 * t92, 0, 0, 0, 0; 0, t90, 0, 0, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end