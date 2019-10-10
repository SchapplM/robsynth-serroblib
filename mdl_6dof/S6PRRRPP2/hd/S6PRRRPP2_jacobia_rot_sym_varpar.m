% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPP2
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
%   Wie in S6PRRRPP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:43
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPP2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:43:00
	% EndTime: 2019-10-09 22:43:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:43:00
	% EndTime: 2019-10-09 22:43:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:43:00
	% EndTime: 2019-10-09 22:43:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:43:00
	% EndTime: 2019-10-09 22:43:00
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-09 22:43:00
	% EndTime: 2019-10-09 22:43:01
	% DurationCPUTime: 0.17s
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
	% StartTime: 2019-10-09 22:43:00
	% EndTime: 2019-10-09 22:43:01
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (714->40), mult. (1993->102), div. (87->9), fcn. (2807->13), ass. (0->57)
	t89 = sin(pkin(6));
	t93 = sin(qJ(3));
	t105 = t89 * t93;
	t91 = cos(pkin(6));
	t94 = sin(qJ(2));
	t103 = t91 * t94;
	t88 = sin(pkin(10));
	t90 = cos(pkin(10));
	t97 = cos(qJ(2));
	t83 = t90 * t103 + t88 * t97;
	t96 = cos(qJ(3));
	t74 = -t90 * t105 + t83 * t96;
	t92 = sin(qJ(4));
	t95 = cos(qJ(4));
	t102 = t91 * t97;
	t98 = t90 * t102 - t88 * t94;
	t66 = t74 * t92 + t98 * t95;
	t104 = t89 * t96;
	t87 = t94 * t104 + t91 * t93;
	t79 = t89 * t97 * t95 + t87 * t92;
	t63 = atan2(-t66, t79);
	t60 = sin(t63);
	t61 = cos(t63);
	t58 = -t60 * t66 + t61 * t79;
	t57 = 0.1e1 / t58 ^ 2;
	t84 = t88 * t102 + t90 * t94;
	t106 = t84 * t95;
	t85 = -t88 * t103 + t90 * t97;
	t76 = t88 * t105 + t85 * t96;
	t69 = t76 * t92 - t106;
	t110 = t57 * t69;
	t70 = t76 * t95 + t84 * t92;
	t65 = 0.1e1 / t70 ^ 2;
	t75 = t88 * t104 - t85 * t93;
	t109 = t65 * t75;
	t78 = 0.1e1 / t79 ^ 2;
	t108 = t66 * t78;
	t107 = t75 ^ 2 * t65;
	t101 = t92 * t96;
	t100 = t92 * t97;
	t99 = -t60 * t79 - t61 * t66;
	t86 = -t94 * t105 + t91 * t96;
	t81 = (t96 * t100 - t94 * t95) * t89;
	t80 = -t89 * t100 + t87 * t95;
	t77 = 0.1e1 / t79;
	t73 = -t90 * t104 - t83 * t93;
	t71 = t98 * t101 - t83 * t95;
	t68 = t74 * t95 - t98 * t92;
	t64 = 0.1e1 / t70;
	t62 = 0.1e1 / (t66 ^ 2 * t78 + 0.1e1);
	t59 = 0.1e1 / (0.1e1 + t107);
	t56 = 0.1e1 / t58;
	t55 = 0.1e1 / (t69 ^ 2 * t57 + 0.1e1);
	t54 = (t86 * t108 - t73 * t77) * t92 * t62;
	t53 = (t81 * t108 - t71 * t77) * t62;
	t52 = (t80 * t108 - t68 * t77) * t62;
	t1 = [0, t53, t54, t52, 0, 0; 0, ((-t84 * t101 - t85 * t95) * t56 - (t99 * t53 - t60 * t71 + t61 * t81) * t110) * t55, (t75 * t92 * t56 - ((-t60 * t73 + t61 * t86) * t92 + t99 * t54) * t110) * t55, (t70 * t56 - (t99 * t52 - t60 * t68 + t61 * t80) * t110) * t55, 0, 0; 0, (t84 * t93 * t64 - (-t96 * t106 + t85 * t92) * t109) * t59, (-t95 * t107 - t64 * t76) * t59, t69 * t59 * t109, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:43:00
	% EndTime: 2019-10-09 22:43:01
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (338->30), mult. (960->79), div. (66->9), fcn. (1365->13), ass. (0->49)
	t81 = sin(pkin(10));
	t83 = cos(pkin(10));
	t90 = cos(qJ(2));
	t84 = cos(pkin(6));
	t87 = sin(qJ(2));
	t93 = t84 * t87;
	t76 = t81 * t90 + t83 * t93;
	t86 = sin(qJ(3));
	t82 = sin(pkin(6));
	t89 = cos(qJ(3));
	t95 = t82 * t89;
	t69 = t76 * t86 + t83 * t95;
	t96 = t82 * t86;
	t79 = t84 * t89 - t87 * t96;
	t68 = atan2(t69, t79);
	t65 = sin(t68);
	t66 = cos(t68);
	t58 = t65 * t69 + t66 * t79;
	t57 = 0.1e1 / t58 ^ 2;
	t78 = -t81 * t93 + t83 * t90;
	t71 = -t78 * t86 + t81 * t95;
	t101 = t57 * t71;
	t72 = t78 * t89 + t81 * t96;
	t92 = t84 * t90;
	t77 = t81 * t92 + t83 * t87;
	t85 = sin(qJ(4));
	t88 = cos(qJ(4));
	t64 = t72 * t88 + t77 * t85;
	t62 = 0.1e1 / t64 ^ 2;
	t63 = -t72 * t85 + t77 * t88;
	t100 = t63 ^ 2 * t62;
	t99 = t62 * t63;
	t74 = 0.1e1 / t79 ^ 2;
	t98 = t69 * t74;
	t97 = t77 * t89;
	t94 = t82 * t90;
	t91 = -t65 * t79 + t66 * t69;
	t80 = -t84 * t86 - t87 * t95;
	t75 = -t81 * t87 + t83 * t92;
	t73 = 0.1e1 / t79;
	t70 = t76 * t89 - t83 * t96;
	t67 = 0.1e1 / (t69 ^ 2 * t74 + 0.1e1);
	t61 = 0.1e1 / t64;
	t59 = 0.1e1 / (0.1e1 + t100);
	t56 = 0.1e1 / t58;
	t55 = 0.1e1 / (t57 * t71 ^ 2 + 0.1e1);
	t54 = (t73 * t75 + t94 * t98) * t86 * t67;
	t53 = (t70 * t73 - t80 * t98) * t67;
	t1 = [0, t54, t53, 0, 0, 0; 0, (t77 * t86 * t56 - ((t65 * t75 - t66 * t94) * t86 + t91 * t54) * t101) * t55, (-t72 * t56 - (t91 * t53 + t65 * t70 + t66 * t80) * t101) * t55, 0, 0, 0; 0, ((t78 * t88 + t85 * t97) * t61 - (t78 * t85 - t88 * t97) * t99) * t59, (-t61 * t85 - t88 * t99) * t71 * t59, (-t61 * t64 - t100) * t59, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end