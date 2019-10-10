% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPP3
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
%   Wie in S6PRRRPP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (101->18), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->29)
	t41 = sin(pkin(10));
	t42 = sin(pkin(6));
	t53 = t41 * t42;
	t48 = cos(qJ(2));
	t52 = t42 * t48;
	t44 = cos(pkin(6));
	t46 = sin(qJ(2));
	t51 = t44 * t46;
	t50 = t44 * t48;
	t43 = cos(pkin(10));
	t37 = -t41 * t51 + t43 * t48;
	t45 = sin(qJ(3));
	t47 = cos(qJ(3));
	t28 = t37 * t47 + t45 * t53;
	t26 = 0.1e1 / t28 ^ 2;
	t27 = t37 * t45 - t47 * t53;
	t49 = t27 ^ 2 * t26 + 0.1e1;
	t40 = 0.1e1 / t48 ^ 2;
	t36 = t41 * t50 + t43 * t46;
	t35 = t41 * t48 + t43 * t51;
	t33 = t41 * t46 - t43 * t50;
	t31 = atan2(-t33, -t52);
	t30 = cos(t31);
	t29 = sin(t31);
	t25 = 0.1e1 / t49;
	t24 = -t29 * t33 - t30 * t52;
	t23 = 0.1e1 / t24 ^ 2;
	t21 = (t35 / t48 + t46 * t33 * t40) / t42 / (0.1e1 + t33 ^ 2 / t42 ^ 2 * t40);
	t1 = [0, t21, 0, 0, 0, 0; 0, (t37 / t24 - (t30 * t42 * t46 - t29 * t35 + (t29 * t52 - t30 * t33) * t21) * t36 * t23) / (t36 ^ 2 * t23 + 0.1e1), 0, 0, 0, 0; 0, (-t45 / t28 + t47 * t27 * t26) * t36 * t25, t49 * t25, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (334->30), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->49)
	t62 = sin(pkin(10));
	t64 = cos(pkin(10));
	t71 = cos(qJ(2));
	t65 = cos(pkin(6));
	t68 = sin(qJ(2));
	t75 = t65 * t68;
	t56 = t62 * t71 + t64 * t75;
	t67 = sin(qJ(3));
	t63 = sin(pkin(6));
	t70 = cos(qJ(3));
	t77 = t63 * t70;
	t48 = t56 * t67 + t64 * t77;
	t78 = t63 * t67;
	t59 = -t65 * t70 + t68 * t78;
	t47 = atan2(-t48, t59);
	t44 = sin(t47);
	t45 = cos(t47);
	t38 = -t44 * t48 + t45 * t59;
	t37 = 0.1e1 / t38 ^ 2;
	t58 = -t62 * t75 + t64 * t71;
	t51 = t58 * t67 - t62 * t77;
	t82 = t37 * t51;
	t52 = t58 * t70 + t62 * t78;
	t74 = t65 * t71;
	t57 = t62 * t74 + t64 * t68;
	t66 = sin(qJ(4));
	t69 = cos(qJ(4));
	t43 = t52 * t69 + t57 * t66;
	t41 = 0.1e1 / t43 ^ 2;
	t42 = t52 * t66 - t57 * t69;
	t81 = t41 * t42;
	t54 = 0.1e1 / t59 ^ 2;
	t80 = t48 * t54;
	t79 = t57 * t70;
	t76 = t63 * t71;
	t73 = t42 ^ 2 * t41 + 0.1e1;
	t72 = -t44 * t59 - t45 * t48;
	t60 = t65 * t67 + t68 * t77;
	t55 = -t62 * t68 + t64 * t74;
	t53 = 0.1e1 / t59;
	t50 = t56 * t70 - t64 * t78;
	t46 = 0.1e1 / (t48 ^ 2 * t54 + 0.1e1);
	t40 = 0.1e1 / t43;
	t39 = 0.1e1 / t73;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / (t51 ^ 2 * t37 + 0.1e1);
	t34 = (-t53 * t55 + t76 * t80) * t67 * t46;
	t33 = (-t50 * t53 + t60 * t80) * t46;
	t1 = [0, t34, t33, 0, 0, 0; 0, (-t57 * t67 * t36 - ((-t44 * t55 + t45 * t76) * t67 + t72 * t34) * t82) * t35, (t52 * t36 - (t72 * t33 - t44 * t50 + t45 * t60) * t82) * t35, 0, 0, 0; 0, ((-t58 * t69 - t66 * t79) * t40 - (t58 * t66 - t69 * t79) * t81) * t39, (-t40 * t66 + t69 * t81) * t51 * t39, t73 * t39, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (707->41), mult. (1976->104), div. (87->9), fcn. (2786->13), ass. (0->56)
	t87 = sin(pkin(6));
	t91 = sin(qJ(3));
	t103 = t87 * t91;
	t89 = cos(pkin(6));
	t92 = sin(qJ(2));
	t101 = t89 * t92;
	t86 = sin(pkin(10));
	t88 = cos(pkin(10));
	t95 = cos(qJ(2));
	t80 = t88 * t101 + t86 * t95;
	t94 = cos(qJ(3));
	t71 = -t88 * t103 + t80 * t94;
	t90 = sin(qJ(4));
	t93 = cos(qJ(4));
	t100 = t89 * t95;
	t96 = t88 * t100 - t86 * t92;
	t62 = t71 * t90 + t96 * t93;
	t102 = t87 * t94;
	t84 = t92 * t102 + t89 * t91;
	t76 = t87 * t95 * t93 + t84 * t90;
	t61 = atan2(-t62, t76);
	t58 = sin(t61);
	t59 = cos(t61);
	t56 = -t58 * t62 + t59 * t76;
	t55 = 0.1e1 / t56 ^ 2;
	t81 = t86 * t100 + t88 * t92;
	t104 = t81 * t93;
	t82 = -t86 * t101 + t88 * t95;
	t73 = t86 * t103 + t82 * t94;
	t65 = t73 * t90 - t104;
	t107 = t55 * t65;
	t75 = 0.1e1 / t76 ^ 2;
	t106 = t62 * t75;
	t66 = t73 * t93 + t81 * t90;
	t72 = -t86 * t102 + t82 * t91;
	t69 = 0.1e1 / t72 ^ 2;
	t105 = t66 * t69;
	t99 = t90 * t94;
	t98 = t90 * t95;
	t97 = -t58 * t76 - t59 * t62;
	t83 = -t92 * t103 + t89 * t94;
	t78 = (-t92 * t93 + t94 * t98) * t87;
	t77 = t84 * t93 - t87 * t98;
	t74 = 0.1e1 / t76;
	t70 = -t88 * t102 - t80 * t91;
	t68 = 0.1e1 / t72;
	t67 = -t80 * t93 + t96 * t99;
	t64 = t71 * t93 - t96 * t90;
	t60 = 0.1e1 / (t62 ^ 2 * t75 + 0.1e1);
	t57 = 0.1e1 / (t66 ^ 2 * t69 + 0.1e1);
	t54 = 0.1e1 / t56;
	t53 = 0.1e1 / (t65 ^ 2 * t55 + 0.1e1);
	t52 = (t83 * t106 - t70 * t74) * t90 * t60;
	t51 = (t78 * t106 - t67 * t74) * t60;
	t50 = (t77 * t106 - t64 * t74) * t60;
	t1 = [0, t51, t52, t50, 0, 0; 0, ((-t81 * t99 - t82 * t93) * t54 - (t97 * t51 - t58 * t67 + t59 * t78) * t107) * t53, (-t72 * t90 * t54 - ((-t58 * t70 + t59 * t83) * t90 + t97 * t52) * t107) * t53, (t66 * t54 - (t97 * t50 - t58 * t64 + t59 * t77) * t107) * t53, 0, 0; 0, ((-t94 * t104 + t82 * t90) * t68 + t81 * t91 * t105) * t57, (-t68 * t72 * t93 - t73 * t105) * t57, -t65 * t68 * t57, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:44:50
	% EndTime: 2019-10-09 22:44:50
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (707->41), mult. (1976->104), div. (87->9), fcn. (2786->13), ass. (0->56)
	t89 = sin(pkin(6));
	t93 = sin(qJ(3));
	t104 = t89 * t93;
	t91 = cos(pkin(6));
	t94 = sin(qJ(2));
	t102 = t91 * t94;
	t88 = sin(pkin(10));
	t90 = cos(pkin(10));
	t97 = cos(qJ(2));
	t82 = t90 * t102 + t88 * t97;
	t96 = cos(qJ(3));
	t73 = -t90 * t104 + t82 * t96;
	t101 = t91 * t97;
	t81 = -t90 * t101 + t88 * t94;
	t92 = sin(qJ(4));
	t95 = cos(qJ(4));
	t65 = t73 * t95 + t81 * t92;
	t103 = t89 * t96;
	t86 = t94 * t103 + t91 * t93;
	t79 = -t89 * t92 * t97 + t86 * t95;
	t63 = atan2(-t65, t79);
	t60 = sin(t63);
	t61 = cos(t63);
	t58 = -t60 * t65 + t61 * t79;
	t57 = 0.1e1 / t58 ^ 2;
	t83 = t88 * t101 + t90 * t94;
	t105 = t83 * t92;
	t84 = -t88 * t102 + t90 * t97;
	t75 = t88 * t104 + t84 * t96;
	t68 = t75 * t95 + t105;
	t108 = t57 * t68;
	t77 = 0.1e1 / t79 ^ 2;
	t107 = t65 * t77;
	t67 = -t75 * t92 + t83 * t95;
	t74 = -t88 * t103 + t84 * t93;
	t71 = 0.1e1 / t74 ^ 2;
	t106 = t67 * t71;
	t100 = t95 * t96;
	t99 = t95 * t97;
	t98 = -t60 * t79 - t61 * t65;
	t85 = -t94 * t104 + t91 * t96;
	t80 = (t92 * t94 + t96 * t99) * t89;
	t78 = -t86 * t92 - t89 * t99;
	t76 = 0.1e1 / t79;
	t72 = -t90 * t103 - t82 * t93;
	t70 = 0.1e1 / t74;
	t69 = -t81 * t100 + t82 * t92;
	t64 = t73 * t92 - t81 * t95;
	t62 = 0.1e1 / (t65 ^ 2 * t77 + 0.1e1);
	t59 = 0.1e1 / (t67 ^ 2 * t71 + 0.1e1);
	t56 = 0.1e1 / t58;
	t55 = 0.1e1 / (t57 * t68 ^ 2 + 0.1e1);
	t54 = (t85 * t107 - t72 * t76) * t95 * t62;
	t53 = (t80 * t107 - t69 * t76) * t62;
	t52 = (t78 * t107 + t64 * t76) * t62;
	t1 = [0, t53, t54, t52, 0, 0; 0, ((-t83 * t100 + t84 * t92) * t56 - (t98 * t53 - t60 * t69 + t61 * t80) * t108) * t55, (-t74 * t95 * t56 - ((-t60 * t72 + t61 * t85) * t95 + t98 * t54) * t108) * t55, (t67 * t56 - (t98 * t52 + t60 * t64 + t61 * t78) * t108) * t55, 0, 0; 0, ((t96 * t105 + t84 * t95) * t70 + t83 * t93 * t106) * t59, (t70 * t74 * t92 - t75 * t106) * t59, -t68 * t70 * t59, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end