% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRP6
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
%   Wie in S6PRPRRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRP6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (68->16), mult. (173->36), div. (41->11), fcn. (271->9), ass. (0->23)
	t34 = sin(pkin(6));
	t38 = cos(qJ(2));
	t41 = t34 * t38;
	t36 = cos(pkin(6));
	t37 = sin(qJ(2));
	t40 = t36 * t37;
	t39 = t36 * t38;
	t35 = cos(pkin(10));
	t33 = sin(pkin(10));
	t32 = 0.1e1 / t38 ^ 2;
	t31 = 0.1e1 / t34 ^ 2;
	t30 = 0.1e1 / t34;
	t28 = -t33 * t40 + t35 * t38;
	t27 = t33 * t39 + t35 * t37;
	t26 = t33 * t38 + t35 * t40;
	t24 = t33 * t37 - t35 * t39;
	t22 = atan2(-t24, -t41);
	t21 = cos(t22);
	t20 = sin(t22);
	t19 = -t20 * t24 - t21 * t41;
	t18 = 0.1e1 / t19 ^ 2;
	t16 = (t26 / t38 + t37 * t24 * t32) * t30 / (t24 ^ 2 * t31 * t32 + 0.1e1);
	t1 = [0, t16, 0, 0, 0, 0; 0, (t28 / t19 - (t21 * t34 * t37 - t20 * t26 + (t20 * t41 - t21 * t24) * t16) * t27 * t18) / (t27 ^ 2 * t18 + 0.1e1), 0, 0, 0, 0; 0, -t27 / t33 * t30 / (0.1e1 + t28 ^ 2 / t33 ^ 2 * t31), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (89->17), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->29)
	t41 = sin(pkin(10));
	t42 = sin(pkin(6));
	t53 = t41 * t42;
	t46 = sin(qJ(2));
	t52 = t42 * t46;
	t44 = cos(pkin(6));
	t51 = t44 * t46;
	t48 = cos(qJ(2));
	t50 = t44 * t48;
	t43 = cos(pkin(10));
	t37 = t41 * t50 + t43 * t46;
	t45 = sin(qJ(4));
	t47 = cos(qJ(4));
	t29 = t37 * t45 + t47 * t53;
	t27 = 0.1e1 / t29 ^ 2;
	t28 = -t37 * t47 + t45 * t53;
	t49 = t28 ^ 2 * t27 + 0.1e1;
	t40 = 0.1e1 / t46 ^ 2;
	t38 = -t41 * t51 + t43 * t48;
	t35 = t41 * t48 + t43 * t51;
	t34 = t41 * t46 - t43 * t50;
	t33 = atan2(-t35, t52);
	t31 = cos(t33);
	t30 = sin(t33);
	t26 = 0.1e1 / t49;
	t25 = -t30 * t35 + t31 * t52;
	t24 = 0.1e1 / t25 ^ 2;
	t22 = (t34 / t46 + t48 * t35 * t40) / t42 / (0.1e1 + t35 ^ 2 / t42 ^ 2 * t40);
	t1 = [0, t22, 0, 0, 0, 0; 0, (-t37 / t25 - (t31 * t42 * t48 + t30 * t34 + (-t30 * t52 - t31 * t35) * t22) * t38 * t24) / (t38 ^ 2 * t24 + 0.1e1), 0, 0, 0, 0; 0, (-t47 / t29 - t45 * t28 * t27) * t38 * t26, 0, t49 * t26, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (334->29), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->51)
	t61 = sin(pkin(10));
	t63 = cos(pkin(10));
	t67 = sin(qJ(2));
	t64 = cos(pkin(6));
	t70 = cos(qJ(2));
	t73 = t64 * t70;
	t55 = t61 * t67 - t63 * t73;
	t69 = cos(qJ(4));
	t62 = sin(pkin(6));
	t66 = sin(qJ(4));
	t78 = t62 * t66;
	t50 = t55 * t69 + t63 * t78;
	t75 = t62 * t70;
	t59 = t64 * t66 + t69 * t75;
	t47 = atan2(t50, t59);
	t44 = sin(t47);
	t45 = cos(t47);
	t38 = t44 * t50 + t45 * t59;
	t37 = 0.1e1 / t38 ^ 2;
	t57 = t61 * t73 + t63 * t67;
	t48 = -t57 * t69 + t61 * t78;
	t83 = t37 * t48;
	t76 = t62 * t69;
	t49 = t57 * t66 + t61 * t76;
	t68 = cos(qJ(5));
	t74 = t64 * t67;
	t58 = -t61 * t74 + t63 * t70;
	t65 = sin(qJ(5));
	t80 = t58 * t65;
	t43 = t49 * t68 + t80;
	t41 = 0.1e1 / t43 ^ 2;
	t79 = t58 * t68;
	t42 = t49 * t65 - t79;
	t82 = t41 * t42;
	t54 = 0.1e1 / t59 ^ 2;
	t81 = t50 * t54;
	t77 = t62 * t67;
	t72 = t42 ^ 2 * t41 + 0.1e1;
	t71 = -t44 * t59 + t45 * t50;
	t60 = t64 * t69 - t66 * t75;
	t56 = t61 * t70 + t63 * t74;
	t53 = 0.1e1 / t59;
	t51 = -t55 * t66 + t63 * t76;
	t46 = 0.1e1 / (t50 ^ 2 * t54 + 0.1e1);
	t40 = 0.1e1 / t43;
	t39 = 0.1e1 / t72;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / (t48 ^ 2 * t37 + 0.1e1);
	t34 = (t53 * t56 + t77 * t81) * t69 * t46;
	t33 = (t51 * t53 - t60 * t81) * t46;
	t1 = [0, t34, 0, t33, 0, 0; 0, (-t58 * t69 * t36 - ((t44 * t56 - t45 * t77) * t69 + t71 * t34) * t83) * t35, 0, (t49 * t36 - (t71 * t33 + t44 * t51 + t45 * t60) * t83) * t35, 0, 0; 0, ((t57 * t68 + t66 * t80) * t40 - (-t57 * t65 + t66 * t79) * t82) * t39, 0, (-t40 * t65 + t68 * t82) * t48 * t39, t72 * t39, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:51:56
	% EndTime: 2019-10-09 21:51:56
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (714->40), mult. (1993->102), div. (87->9), fcn. (2807->13), ass. (0->59)
	t92 = cos(pkin(6));
	t95 = sin(qJ(2));
	t104 = t92 * t95;
	t89 = sin(pkin(10));
	t91 = cos(pkin(10));
	t98 = cos(qJ(2));
	t100 = t91 * t104 + t89 * t98;
	t90 = sin(pkin(6));
	t97 = cos(qJ(4));
	t106 = t90 * t97;
	t103 = t92 * t98;
	t84 = -t91 * t103 + t89 * t95;
	t94 = sin(qJ(4));
	t77 = -t91 * t106 + t84 * t94;
	t93 = sin(qJ(5));
	t96 = cos(qJ(5));
	t69 = -t100 * t96 + t77 * t93;
	t105 = t90 * t98;
	t88 = -t94 * t105 + t92 * t97;
	t80 = -t90 * t95 * t96 + t88 * t93;
	t64 = atan2(-t69, t80);
	t61 = sin(t64);
	t62 = cos(t64);
	t59 = -t61 * t69 + t62 * t80;
	t58 = 0.1e1 / t59 ^ 2;
	t86 = -t89 * t104 + t91 * t98;
	t108 = t86 * t96;
	t85 = t89 * t103 + t91 * t95;
	t75 = t89 * t106 + t85 * t94;
	t67 = t75 * t93 - t108;
	t113 = t58 * t67;
	t109 = t86 * t93;
	t68 = t75 * t96 + t109;
	t66 = 0.1e1 / t68 ^ 2;
	t107 = t90 * t94;
	t74 = -t89 * t107 + t85 * t97;
	t112 = t66 * t74;
	t79 = 0.1e1 / t80 ^ 2;
	t111 = t69 * t79;
	t110 = t74 ^ 2 * t66;
	t102 = t93 * t95;
	t101 = -t61 * t80 - t62 * t69;
	t99 = t100 * t93;
	t87 = -t97 * t105 - t92 * t94;
	t82 = (t94 * t102 - t96 * t98) * t90;
	t81 = t90 * t102 + t88 * t96;
	t78 = 0.1e1 / t80;
	t76 = t91 * t107 + t84 * t97;
	t72 = t84 * t96 + t94 * t99;
	t71 = t77 * t96 + t99;
	t65 = 0.1e1 / t68;
	t63 = 0.1e1 / (t69 ^ 2 * t79 + 0.1e1);
	t60 = 0.1e1 / (0.1e1 + t110);
	t57 = 0.1e1 / t59;
	t56 = 0.1e1 / (t67 ^ 2 * t58 + 0.1e1);
	t55 = (t87 * t111 - t76 * t78) * t93 * t63;
	t54 = (t82 * t111 - t72 * t78) * t63;
	t53 = (t81 * t111 - t71 * t78) * t63;
	t1 = [0, t54, 0, t55, t53, 0; 0, ((t94 * t109 + t85 * t96) * t57 - (t101 * t54 - t61 * t72 + t62 * t82) * t113) * t56, 0, (t74 * t93 * t57 - ((-t61 * t76 + t62 * t87) * t93 + t101 * t55) * t113) * t56, (t68 * t57 - (t101 * t53 - t61 * t71 + t62 * t81) * t113) * t56, 0; 0, (t86 * t97 * t65 - (t94 * t108 - t85 * t93) * t112) * t60, 0, (-t96 * t110 - t65 * t75) * t60, t67 * t60 * t112, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end