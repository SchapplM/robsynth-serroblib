% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR6
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
%   Wie in S6RRRPPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:26
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (35->13), mult. (116->34), div. (25->9), fcn. (175->9), ass. (0->25)
	t30 = cos(pkin(6));
	t29 = sin(pkin(6));
	t34 = cos(qJ(1));
	t38 = t34 * t29;
	t26 = atan2(t38, t30);
	t23 = sin(t26);
	t24 = cos(t26);
	t18 = t23 * t38 + t24 * t30;
	t32 = sin(qJ(1));
	t42 = 0.1e1 / t18 ^ 2 * t32 ^ 2;
	t27 = t29 ^ 2;
	t25 = 0.1e1 / (0.1e1 + t34 ^ 2 * t27 / t30 ^ 2);
	t41 = t25 / t30;
	t31 = sin(qJ(2));
	t40 = t32 * t31;
	t33 = cos(qJ(2));
	t39 = t32 * t33;
	t37 = t34 * t31;
	t36 = t34 * t33;
	t22 = -t30 * t40 + t36;
	t20 = 0.1e1 / t22 ^ 2;
	t21 = t30 * t39 + t37;
	t35 = t21 ^ 2 * t20 + 0.1e1;
	t19 = 0.1e1 / t35;
	t1 = [-t32 * t29 * t41, 0, 0, 0, 0, 0; (0.1e1 / t18 * t38 - (-t24 * t27 * t34 * t41 + (t25 - 0.1e1) * t29 * t23) * t29 * t42) / (t27 * t42 + 0.1e1), 0, 0, 0, 0, 0; ((t30 * t36 - t40) / t22 - (-t30 * t37 - t39) * t21 * t20) * t19, t35 * t19, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (173->24), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->40)
	t53 = sin(qJ(2));
	t56 = cos(qJ(2));
	t57 = cos(qJ(1));
	t54 = sin(qJ(1));
	t61 = cos(pkin(6));
	t59 = t54 * t61;
	t46 = -t53 * t59 + t57 * t56;
	t52 = sin(qJ(3));
	t55 = cos(qJ(3));
	t51 = sin(pkin(6));
	t64 = t51 * t54;
	t37 = t46 * t55 + t52 * t64;
	t35 = 0.1e1 / t37 ^ 2;
	t36 = t46 * t52 - t55 * t64;
	t68 = t35 * t36;
	t58 = t57 * t61;
	t42 = t54 * t53 - t56 * t58;
	t63 = t51 * t56;
	t40 = atan2(-t42, -t63);
	t39 = cos(t40);
	t67 = t39 * t42;
	t38 = sin(t40);
	t32 = -t38 * t42 - t39 * t63;
	t31 = 0.1e1 / t32 ^ 2;
	t45 = t57 * t53 + t56 * t59;
	t66 = t45 ^ 2 * t31;
	t48 = 0.1e1 / t51;
	t49 = 0.1e1 / t56;
	t65 = t48 * t49;
	t62 = t51 * t57;
	t60 = t36 ^ 2 * t35 + 0.1e1;
	t50 = 0.1e1 / t56 ^ 2;
	t44 = t53 * t58 + t54 * t56;
	t41 = 0.1e1 / (0.1e1 + t42 ^ 2 / t51 ^ 2 * t50);
	t34 = 0.1e1 / t37;
	t33 = 0.1e1 / t60;
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / (0.1e1 + t66);
	t28 = (t42 * t50 * t53 + t44 * t49) * t48 * t41;
	t1 = [t45 * t41 * t65, t28, 0, 0, 0, 0; (-t42 * t30 - (-t38 + (-t65 * t67 + t38) * t41) * t66) * t29, (t46 * t30 - (t39 * t51 * t53 - t38 * t44 + (t38 * t63 - t67) * t28) * t45 * t31) * t29, 0, 0, 0, 0; ((-t44 * t52 - t55 * t62) * t34 - (-t44 * t55 + t52 * t62) * t68) * t33, (-t52 * t34 + t55 * t68) * t45 * t33, t60 * t33, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (221->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
	t59 = sin(qJ(2));
	t61 = cos(qJ(2));
	t62 = cos(qJ(1));
	t60 = sin(qJ(1));
	t66 = cos(pkin(6));
	t64 = t60 * t66;
	t50 = -t59 * t64 + t62 * t61;
	t55 = qJ(3) + pkin(11);
	t52 = sin(t55);
	t53 = cos(t55);
	t58 = sin(pkin(6));
	t69 = t58 * t60;
	t41 = t50 * t53 + t52 * t69;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = t50 * t52 - t53 * t69;
	t73 = t39 * t40;
	t63 = t62 * t66;
	t46 = t60 * t59 - t61 * t63;
	t68 = t58 * t61;
	t44 = atan2(-t46, -t68);
	t43 = cos(t44);
	t72 = t43 * t46;
	t42 = sin(t44);
	t36 = -t42 * t46 - t43 * t68;
	t35 = 0.1e1 / t36 ^ 2;
	t49 = t62 * t59 + t61 * t64;
	t71 = t49 ^ 2 * t35;
	t54 = 0.1e1 / t58;
	t56 = 0.1e1 / t61;
	t70 = t54 * t56;
	t67 = t58 * t62;
	t65 = t40 ^ 2 * t39 + 0.1e1;
	t57 = 0.1e1 / t61 ^ 2;
	t48 = t59 * t63 + t60 * t61;
	t45 = 0.1e1 / (0.1e1 + t46 ^ 2 / t58 ^ 2 * t57);
	t38 = 0.1e1 / t41;
	t37 = 0.1e1 / t65;
	t34 = 0.1e1 / t36;
	t33 = 0.1e1 / (0.1e1 + t71);
	t32 = (t46 * t57 * t59 + t48 * t56) * t54 * t45;
	t1 = [t49 * t45 * t70, t32, 0, 0, 0, 0; (-t46 * t34 - (-t42 + (-t70 * t72 + t42) * t45) * t71) * t33, (t50 * t34 - (t43 * t58 * t59 - t42 * t48 + (t42 * t68 - t72) * t32) * t49 * t35) * t33, 0, 0, 0, 0; ((-t48 * t52 - t53 * t67) * t38 - (-t48 * t53 + t52 * t67) * t73) * t37, (-t52 * t38 + t53 * t73) * t49 * t37, t65 * t37, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (727->32), mult. (1022->76), div. (77->9), fcn. (1479->11), ass. (0->49)
	t74 = cos(pkin(6));
	t75 = sin(qJ(2));
	t78 = cos(qJ(1));
	t81 = t78 * t75;
	t76 = sin(qJ(1));
	t77 = cos(qJ(2));
	t82 = t76 * t77;
	t65 = t74 * t81 + t82;
	t72 = qJ(3) + pkin(11);
	t70 = sin(t72);
	t71 = cos(t72);
	t73 = sin(pkin(6));
	t84 = t73 * t78;
	t53 = t65 * t70 + t71 * t84;
	t87 = t73 * t75;
	t60 = t70 * t87 - t74 * t71;
	t51 = atan2(-t53, t60);
	t48 = sin(t51);
	t49 = cos(t51);
	t47 = -t48 * t53 + t49 * t60;
	t46 = 0.1e1 / t47 ^ 2;
	t80 = t78 * t77;
	t83 = t76 * t75;
	t67 = -t74 * t83 + t80;
	t86 = t73 * t76;
	t56 = t67 * t70 - t71 * t86;
	t92 = t46 * t56;
	t91 = t49 * t53;
	t59 = 0.1e1 / t60 ^ 2;
	t90 = t53 * t59;
	t89 = t56 ^ 2 * t46;
	t57 = t67 * t71 + t70 * t86;
	t66 = t74 * t82 + t81;
	t63 = 0.1e1 / t66 ^ 2;
	t88 = t57 * t63;
	t85 = t73 * t77;
	t55 = t65 * t71 - t70 * t84;
	t79 = -t48 * t60 - t91;
	t64 = t74 * t80 - t83;
	t62 = 0.1e1 / t66;
	t61 = t74 * t70 + t71 * t87;
	t58 = 0.1e1 / t60;
	t52 = 0.1e1 / (t57 ^ 2 * t63 + 0.1e1);
	t50 = 0.1e1 / (t53 ^ 2 * t59 + 0.1e1);
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / (0.1e1 + t89);
	t43 = (-t58 * t64 + t85 * t90) * t70 * t50;
	t42 = (-t55 * t58 + t61 * t90) * t50;
	t1 = [-t56 * t58 * t50, t43, t42, 0, 0, 0; (-t53 * t45 - (-t48 + (t58 * t91 + t48) * t50) * t89) * t44, (-t66 * t70 * t45 - ((-t48 * t64 + t49 * t85) * t70 + t79 * t43) * t92) * t44, (t57 * t45 - (t42 * t79 - t48 * t55 + t49 * t61) * t92) * t44, 0, 0, 0; (-t55 * t62 - t64 * t88) * t52, (-t62 * t66 * t71 - t67 * t88) * t52, -t56 * t62 * t52, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:25:56
	% EndTime: 2019-10-10 11:25:56
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (878->38), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
	t81 = cos(pkin(6));
	t84 = sin(qJ(1));
	t86 = cos(qJ(2));
	t91 = t84 * t86;
	t83 = sin(qJ(2));
	t87 = cos(qJ(1));
	t93 = t83 * t87;
	t73 = t81 * t93 + t91;
	t79 = qJ(3) + pkin(11);
	t77 = sin(t79);
	t78 = cos(t79);
	t80 = sin(pkin(6));
	t94 = t80 * t87;
	t63 = t73 * t78 - t77 * t94;
	t97 = t80 * t83;
	t71 = t77 * t81 + t78 * t97;
	t61 = atan2(-t63, t71);
	t54 = sin(t61);
	t55 = cos(t61);
	t52 = -t54 * t63 + t55 * t71;
	t51 = 0.1e1 / t52 ^ 2;
	t90 = t86 * t87;
	t92 = t84 * t83;
	t75 = -t81 * t92 + t90;
	t96 = t80 * t84;
	t67 = t75 * t78 + t77 * t96;
	t104 = t51 * t67;
	t103 = t51 * t67 ^ 2;
	t102 = t55 * t63;
	t66 = t75 * t77 - t78 * t96;
	t82 = sin(qJ(6));
	t74 = t81 * t91 + t93;
	t85 = cos(qJ(6));
	t98 = t74 * t85;
	t60 = t66 * t82 + t98;
	t57 = 0.1e1 / t60 ^ 2;
	t99 = t74 * t82;
	t59 = -t66 * t85 + t99;
	t101 = t57 * t59;
	t69 = 0.1e1 / t71 ^ 2;
	t100 = t63 * t69;
	t95 = t80 * t86;
	t89 = t57 * t59 ^ 2 + 0.1e1;
	t88 = -t54 * t71 - t102;
	t62 = t73 * t77 + t78 * t94;
	t72 = t81 * t90 - t92;
	t70 = -t77 * t97 + t78 * t81;
	t68 = 0.1e1 / t71;
	t58 = 0.1e1 / (t63 ^ 2 * t69 + 0.1e1);
	t56 = 0.1e1 / t60;
	t53 = 0.1e1 / t89;
	t50 = 0.1e1 / t52;
	t49 = 0.1e1 / (0.1e1 + t103);
	t48 = (t95 * t100 - t68 * t72) * t78 * t58;
	t47 = (t70 * t100 + t62 * t68) * t58;
	t1 = [-t67 * t68 * t58, t48, t47, 0, 0, 0; (-t63 * t50 - (-t54 + (t68 * t102 + t54) * t58) * t103) * t49, (-t74 * t78 * t50 - ((-t54 * t72 + t55 * t95) * t78 + t88 * t48) * t104) * t49, (-t66 * t50 - (t88 * t47 + t54 * t62 + t55 * t70) * t104) * t49, 0, 0, 0; ((t62 * t85 + t72 * t82) * t56 - (-t62 * t82 + t72 * t85) * t101) * t53, ((t75 * t82 + t77 * t98) * t56 - (t75 * t85 - t77 * t99) * t101) * t53, (-t82 * t101 - t85 * t56) * t67 * t53, 0, 0, t89 * t53;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end