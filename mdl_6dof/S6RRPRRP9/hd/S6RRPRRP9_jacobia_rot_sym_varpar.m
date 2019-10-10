% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP9
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
%   Wie in S6RRPRRP9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:42
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP9_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP9_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:53
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:53
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:53
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:53
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (156->24), mult. (403->63), div. (67->11), fcn. (610->11), ass. (0->39)
	t46 = sin(qJ(2));
	t48 = cos(qJ(2));
	t49 = cos(qJ(1));
	t47 = sin(qJ(1));
	t52 = cos(pkin(6));
	t51 = t47 * t52;
	t38 = -t46 * t51 + t49 * t48;
	t43 = sin(pkin(11));
	t45 = cos(pkin(11));
	t44 = sin(pkin(6));
	t55 = t44 * t47;
	t29 = t38 * t45 + t43 * t55;
	t27 = 0.1e1 / t29 ^ 2;
	t28 = t38 * t43 - t45 * t55;
	t59 = t27 * t28;
	t50 = t49 * t52;
	t34 = t47 * t46 - t48 * t50;
	t54 = t44 * t48;
	t32 = atan2(-t34, -t54);
	t31 = cos(t32);
	t58 = t31 * t34;
	t30 = sin(t32);
	t24 = -t30 * t34 - t31 * t54;
	t23 = 0.1e1 / t24 ^ 2;
	t37 = t49 * t46 + t48 * t51;
	t57 = t37 ^ 2 * t23;
	t40 = 0.1e1 / t44;
	t41 = 0.1e1 / t48;
	t56 = t40 * t41;
	t53 = t44 * t49;
	t42 = 0.1e1 / t48 ^ 2;
	t36 = t46 * t50 + t47 * t48;
	t33 = 0.1e1 / (0.1e1 + t34 ^ 2 / t44 ^ 2 * t42);
	t26 = 0.1e1 / t29;
	t25 = 0.1e1 / (t28 ^ 2 * t27 + 0.1e1);
	t22 = 0.1e1 / t24;
	t21 = 0.1e1 / (0.1e1 + t57);
	t20 = (t34 * t42 * t46 + t36 * t41) * t40 * t33;
	t1 = [t37 * t33 * t56, t20, 0, 0, 0, 0; (-t34 * t22 - (-t30 + (-t56 * t58 + t30) * t33) * t57) * t21, (t38 * t22 - (t31 * t44 * t46 - t30 * t36 + (t30 * t54 - t58) * t20) * t37 * t23) * t21, 0, 0, 0, 0; ((-t36 * t43 - t45 * t53) * t26 - (-t36 * t45 + t43 * t53) * t59) * t25, (-t43 * t26 + t45 * t59) * t37 * t25, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:53
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (221->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
	t58 = sin(qJ(2));
	t60 = cos(qJ(2));
	t61 = cos(qJ(1));
	t59 = sin(qJ(1));
	t65 = cos(pkin(6));
	t63 = t59 * t65;
	t49 = -t58 * t63 + t61 * t60;
	t54 = pkin(11) + qJ(4);
	t51 = sin(t54);
	t52 = cos(t54);
	t57 = sin(pkin(6));
	t68 = t57 * t59;
	t40 = t49 * t52 + t51 * t68;
	t38 = 0.1e1 / t40 ^ 2;
	t39 = t49 * t51 - t52 * t68;
	t72 = t38 * t39;
	t62 = t61 * t65;
	t45 = t59 * t58 - t60 * t62;
	t67 = t57 * t60;
	t43 = atan2(-t45, -t67);
	t42 = cos(t43);
	t71 = t42 * t45;
	t41 = sin(t43);
	t35 = -t41 * t45 - t42 * t67;
	t34 = 0.1e1 / t35 ^ 2;
	t48 = t61 * t58 + t60 * t63;
	t70 = t48 ^ 2 * t34;
	t53 = 0.1e1 / t57;
	t55 = 0.1e1 / t60;
	t69 = t53 * t55;
	t66 = t57 * t61;
	t64 = t39 ^ 2 * t38 + 0.1e1;
	t56 = 0.1e1 / t60 ^ 2;
	t47 = t58 * t62 + t59 * t60;
	t44 = 0.1e1 / (0.1e1 + t45 ^ 2 / t57 ^ 2 * t56);
	t37 = 0.1e1 / t40;
	t36 = 0.1e1 / t64;
	t33 = 0.1e1 / t35;
	t32 = 0.1e1 / (0.1e1 + t70);
	t31 = (t45 * t56 * t58 + t47 * t55) * t53 * t44;
	t1 = [t48 * t44 * t69, t31, 0, 0, 0, 0; (-t45 * t33 - (-t41 + (-t69 * t71 + t41) * t44) * t70) * t32, (t49 * t33 - (t42 * t57 * t58 - t41 * t47 + (t41 * t67 - t71) * t31) * t48 * t34) * t32, 0, 0, 0, 0; ((-t47 * t51 - t52 * t66) * t37 - (-t47 * t52 + t51 * t66) * t72) * t36, (-t51 * t37 + t52 * t72) * t48 * t36, 0, t64 * t36, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:54
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (878->38), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
	t84 = cos(pkin(6));
	t86 = sin(qJ(2));
	t90 = cos(qJ(1));
	t94 = t90 * t86;
	t87 = sin(qJ(1));
	t89 = cos(qJ(2));
	t95 = t87 * t89;
	t75 = t84 * t94 + t95;
	t82 = pkin(11) + qJ(4);
	t80 = sin(t82);
	t81 = cos(t82);
	t83 = sin(pkin(6));
	t97 = t83 * t90;
	t64 = t75 * t80 + t81 * t97;
	t100 = t83 * t86;
	t72 = t80 * t100 - t84 * t81;
	t63 = atan2(-t64, t72);
	t56 = sin(t63);
	t57 = cos(t63);
	t54 = -t56 * t64 + t57 * t72;
	t53 = 0.1e1 / t54 ^ 2;
	t93 = t90 * t89;
	t96 = t87 * t86;
	t77 = -t84 * t96 + t93;
	t99 = t83 * t87;
	t68 = t77 * t80 - t81 * t99;
	t107 = t53 * t68;
	t106 = t57 * t64;
	t76 = t84 * t95 + t94;
	t85 = sin(qJ(5));
	t102 = t76 * t85;
	t69 = t77 * t81 + t80 * t99;
	t88 = cos(qJ(5));
	t62 = t69 * t88 + t102;
	t59 = 0.1e1 / t62 ^ 2;
	t101 = t76 * t88;
	t61 = t69 * t85 - t101;
	t105 = t59 * t61;
	t71 = 0.1e1 / t72 ^ 2;
	t104 = t64 * t71;
	t103 = t68 ^ 2 * t53;
	t98 = t83 * t89;
	t92 = t61 ^ 2 * t59 + 0.1e1;
	t66 = t75 * t81 - t80 * t97;
	t91 = -t56 * t72 - t106;
	t74 = t84 * t93 - t96;
	t73 = t81 * t100 + t84 * t80;
	t70 = 0.1e1 / t72;
	t60 = 0.1e1 / (t64 ^ 2 * t71 + 0.1e1);
	t58 = 0.1e1 / t62;
	t55 = 0.1e1 / t92;
	t52 = 0.1e1 / t54;
	t51 = 0.1e1 / (0.1e1 + t103);
	t50 = (t98 * t104 - t70 * t74) * t80 * t60;
	t49 = (t73 * t104 - t66 * t70) * t60;
	t1 = [-t68 * t70 * t60, t50, 0, t49, 0, 0; (-t64 * t52 - (-t56 + (t70 * t106 + t56) * t60) * t103) * t51, (-t76 * t80 * t52 - ((-t56 * t74 + t57 * t98) * t80 + t91 * t50) * t107) * t51, 0, (t69 * t52 - (t91 * t49 - t56 * t66 + t57 * t73) * t107) * t51, 0, 0; ((-t66 * t85 - t74 * t88) * t58 - (-t66 * t88 + t74 * t85) * t105) * t55, ((-t81 * t102 - t77 * t88) * t58 - (-t81 * t101 + t77 * t85) * t105) * t55, 0, (t88 * t105 - t85 * t58) * t68 * t55, t92 * t55, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:42:53
	% EndTime: 2019-10-10 10:42:54
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (878->38), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
	t86 = cos(pkin(6));
	t88 = sin(qJ(2));
	t92 = cos(qJ(1));
	t96 = t92 * t88;
	t89 = sin(qJ(1));
	t91 = cos(qJ(2));
	t97 = t89 * t91;
	t77 = t86 * t96 + t97;
	t84 = pkin(11) + qJ(4);
	t82 = sin(t84);
	t83 = cos(t84);
	t85 = sin(pkin(6));
	t99 = t85 * t92;
	t66 = t77 * t82 + t83 * t99;
	t102 = t85 * t88;
	t74 = t82 * t102 - t86 * t83;
	t65 = atan2(-t66, t74);
	t58 = sin(t65);
	t59 = cos(t65);
	t56 = -t58 * t66 + t59 * t74;
	t55 = 0.1e1 / t56 ^ 2;
	t101 = t85 * t89;
	t95 = t92 * t91;
	t98 = t89 * t88;
	t79 = -t86 * t98 + t95;
	t70 = -t83 * t101 + t79 * t82;
	t109 = t55 * t70;
	t108 = t59 * t66;
	t78 = t86 * t97 + t96;
	t87 = sin(qJ(5));
	t104 = t78 * t87;
	t71 = t82 * t101 + t79 * t83;
	t90 = cos(qJ(5));
	t64 = t71 * t90 + t104;
	t61 = 0.1e1 / t64 ^ 2;
	t103 = t78 * t90;
	t63 = t71 * t87 - t103;
	t107 = t61 * t63;
	t73 = 0.1e1 / t74 ^ 2;
	t106 = t66 * t73;
	t105 = t70 ^ 2 * t55;
	t100 = t85 * t91;
	t94 = t63 ^ 2 * t61 + 0.1e1;
	t68 = t77 * t83 - t82 * t99;
	t93 = -t58 * t74 - t108;
	t76 = t86 * t95 - t98;
	t75 = t83 * t102 + t86 * t82;
	t72 = 0.1e1 / t74;
	t62 = 0.1e1 / (t66 ^ 2 * t73 + 0.1e1);
	t60 = 0.1e1 / t64;
	t57 = 0.1e1 / t94;
	t54 = 0.1e1 / t56;
	t53 = 0.1e1 / (0.1e1 + t105);
	t52 = (t100 * t106 - t72 * t76) * t82 * t62;
	t51 = (t75 * t106 - t68 * t72) * t62;
	t1 = [-t70 * t72 * t62, t52, 0, t51, 0, 0; (-t66 * t54 - (-t58 + (t72 * t108 + t58) * t62) * t105) * t53, (-t78 * t82 * t54 - ((t59 * t100 - t58 * t76) * t82 + t93 * t52) * t109) * t53, 0, (t71 * t54 - (t93 * t51 - t58 * t68 + t59 * t75) * t109) * t53, 0, 0; ((-t68 * t87 - t76 * t90) * t60 - (-t68 * t90 + t76 * t87) * t107) * t57, ((-t83 * t104 - t79 * t90) * t60 - (-t83 * t103 + t79 * t87) * t107) * t57, 0, (t90 * t107 - t87 * t60) * t70 * t57, t94 * t57, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end