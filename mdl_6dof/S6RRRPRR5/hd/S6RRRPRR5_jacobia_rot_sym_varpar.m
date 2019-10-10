% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RRRPRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:00
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (297->18), mult. (231->39), div. (67->9), fcn. (349->7), ass. (0->30)
	t55 = cos(qJ(1));
	t52 = t55 ^ 2;
	t53 = qJ(2) + qJ(3);
	t49 = cos(t53);
	t48 = sin(t53);
	t54 = sin(qJ(1));
	t59 = t54 * t48;
	t42 = atan2(-t59, -t49);
	t40 = sin(t42);
	t41 = cos(t42);
	t38 = -t40 * t59 - t41 * t49;
	t37 = 0.1e1 / t38 ^ 2;
	t65 = t37 * t48;
	t64 = t40 * t49;
	t45 = t48 ^ 2;
	t56 = t49 ^ 2;
	t63 = t45 / t56;
	t62 = t48 * t55;
	t57 = t54 ^ 2;
	t61 = 0.1e1 / t57 * t52;
	t43 = 0.1e1 / (t57 * t63 + 0.1e1);
	t60 = t54 * t43;
	t44 = 0.1e1 / (t56 * t61 + 0.1e1);
	t58 = 0.1e1 / t54 * t44 * t62;
	t46 = 0.1e1 / t49;
	t39 = (0.1e1 + t63) * t60;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / (t52 * t45 * t37 + 0.1e1);
	t34 = (t49 * t36 - (-t54 * t64 + t41 * t48 + (-t41 * t59 + t64) * t39) * t65) * t55 * t35;
	t1 = [t46 * t43 * t62, t39, t39, 0, 0, 0; (-t36 * t59 - (-t41 * t45 * t46 * t60 + (t43 - 0.1e1) * t48 * t40) * t52 * t65) * t35, t34, t34, 0, 0, 0; (-0.1e1 - t61) * t49 * t44, -t58, -t58, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (324->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t70 = qJ(2) + qJ(3);
	t68 = sin(t70);
	t69 = cos(t70);
	t72 = sin(qJ(1));
	t80 = t72 * t69;
	t63 = atan2(-t80, t68);
	t61 = sin(t63);
	t62 = cos(t63);
	t54 = -t61 * t80 + t62 * t68;
	t53 = 0.1e1 / t54 ^ 2;
	t74 = cos(qJ(1));
	t86 = t53 * t74 ^ 2;
	t71 = sin(qJ(5));
	t77 = t74 * t71;
	t73 = cos(qJ(5));
	t78 = t72 * t73;
	t60 = t68 * t77 + t78;
	t58 = 0.1e1 / t60 ^ 2;
	t76 = t74 * t73;
	t79 = t72 * t71;
	t59 = -t68 * t76 + t79;
	t85 = t58 * t59;
	t84 = t61 * t68;
	t67 = t69 ^ 2;
	t83 = 0.1e1 / t68 ^ 2 * t67;
	t82 = t69 * t74;
	t64 = 0.1e1 / (t72 ^ 2 * t83 + 0.1e1);
	t81 = t72 * t64;
	t75 = t59 ^ 2 * t58 + 0.1e1;
	t65 = 0.1e1 / t68;
	t57 = 0.1e1 / t60;
	t56 = 0.1e1 / t75;
	t55 = (0.1e1 + t83) * t81;
	t52 = 0.1e1 / t54;
	t51 = 0.1e1 / (t67 * t86 + 0.1e1);
	t50 = (-t57 * t73 - t71 * t85) * t56 * t82;
	t49 = (-t68 * t52 - (t72 * t84 + t62 * t69 + (-t62 * t80 - t84) * t55) * t69 * t53) * t74 * t51;
	t1 = [-t65 * t64 * t82, t55, t55, 0, 0, 0; (-t52 * t80 - (t62 * t65 * t67 * t81 + (t64 - 0.1e1) * t69 * t61) * t69 * t86) * t51, t49, t49, 0, 0, 0; ((t68 * t78 + t77) * t57 - (-t68 * t79 + t76) * t85) * t56, t50, t50, 0, t75 * t56, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:00:36
	% EndTime: 2019-10-10 12:00:36
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (419->22), mult. (359->54), div. (84->9), fcn. (524->9), ass. (0->40)
	t94 = qJ(2) + qJ(3);
	t92 = cos(t94);
	t95 = sin(qJ(1));
	t100 = t95 * t92;
	t90 = sin(t94);
	t84 = atan2(-t100, t90);
	t82 = sin(t84);
	t83 = cos(t84);
	t75 = -t100 * t82 + t83 * t90;
	t74 = 0.1e1 / t75 ^ 2;
	t96 = cos(qJ(1));
	t108 = t74 * t96 ^ 2;
	t93 = qJ(5) + qJ(6);
	t91 = cos(t93);
	t101 = t95 * t91;
	t89 = sin(t93);
	t99 = t96 * t89;
	t81 = t90 * t99 + t101;
	t79 = 0.1e1 / t81 ^ 2;
	t102 = t95 * t89;
	t98 = t96 * t91;
	t80 = -t90 * t98 + t102;
	t107 = t79 * t80;
	t106 = t82 * t90;
	t88 = t92 ^ 2;
	t105 = 0.1e1 / t90 ^ 2 * t88;
	t104 = t92 * t96;
	t85 = 0.1e1 / (t95 ^ 2 * t105 + 0.1e1);
	t103 = t95 * t85;
	t97 = t80 ^ 2 * t79 + 0.1e1;
	t86 = 0.1e1 / t90;
	t78 = 0.1e1 / t81;
	t77 = (0.1e1 + t105) * t103;
	t76 = 0.1e1 / t97;
	t73 = 0.1e1 / t75;
	t72 = 0.1e1 / (t88 * t108 + 0.1e1);
	t71 = t97 * t76;
	t70 = (-t89 * t107 - t78 * t91) * t76 * t104;
	t69 = (-t90 * t73 - (t95 * t106 + t83 * t92 + (-t100 * t83 - t106) * t77) * t92 * t74) * t96 * t72;
	t1 = [-t86 * t85 * t104, t77, t77, 0, 0, 0; (-t73 * t100 - (t83 * t86 * t88 * t103 + (t85 - 0.1e1) * t92 * t82) * t92 * t108) * t72, t69, t69, 0, 0, 0; ((t101 * t90 + t99) * t78 - (-t102 * t90 + t98) * t107) * t76, t70, t70, 0, t71, t71;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end