% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR10
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
%   Wie in S6RRPRPR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:20
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR10_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR10_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
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
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
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
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (727->32), mult. (1022->76), div. (77->9), fcn. (1479->11), ass. (0->49)
	t71 = cos(pkin(6));
	t72 = sin(qJ(2));
	t75 = cos(qJ(1));
	t78 = t75 * t72;
	t73 = sin(qJ(1));
	t74 = cos(qJ(2));
	t79 = t73 * t74;
	t62 = t71 * t78 + t79;
	t69 = pkin(11) + qJ(4);
	t67 = sin(t69);
	t68 = cos(t69);
	t70 = sin(pkin(6));
	t81 = t70 * t75;
	t50 = t62 * t67 + t68 * t81;
	t84 = t70 * t72;
	t57 = t67 * t84 - t71 * t68;
	t48 = atan2(-t50, t57);
	t45 = sin(t48);
	t46 = cos(t48);
	t44 = -t45 * t50 + t46 * t57;
	t43 = 0.1e1 / t44 ^ 2;
	t77 = t75 * t74;
	t80 = t73 * t72;
	t64 = -t71 * t80 + t77;
	t83 = t70 * t73;
	t53 = t64 * t67 - t68 * t83;
	t89 = t43 * t53;
	t88 = t46 * t50;
	t56 = 0.1e1 / t57 ^ 2;
	t87 = t50 * t56;
	t86 = t53 ^ 2 * t43;
	t54 = t64 * t68 + t67 * t83;
	t63 = t71 * t79 + t78;
	t60 = 0.1e1 / t63 ^ 2;
	t85 = t54 * t60;
	t82 = t70 * t74;
	t52 = t62 * t68 - t67 * t81;
	t76 = -t45 * t57 - t88;
	t61 = t71 * t77 - t80;
	t59 = 0.1e1 / t63;
	t58 = t71 * t67 + t68 * t84;
	t55 = 0.1e1 / t57;
	t49 = 0.1e1 / (t54 ^ 2 * t60 + 0.1e1);
	t47 = 0.1e1 / (t50 ^ 2 * t56 + 0.1e1);
	t42 = 0.1e1 / t44;
	t41 = 0.1e1 / (0.1e1 + t86);
	t40 = (-t55 * t61 + t82 * t87) * t67 * t47;
	t39 = (-t52 * t55 + t58 * t87) * t47;
	t1 = [-t53 * t55 * t47, t40, 0, t39, 0, 0; (-t50 * t42 - (-t45 + (t55 * t88 + t45) * t47) * t86) * t41, (-t63 * t67 * t42 - ((-t45 * t61 + t46 * t82) * t67 + t76 * t40) * t89) * t41, 0, (t54 * t42 - (t39 * t76 - t45 * t52 + t46 * t58) * t89) * t41, 0, 0; (-t52 * t59 - t61 * t85) * t49, (-t59 * t63 * t68 - t64 * t85) * t49, 0, -t53 * t59 * t49, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:20:52
	% EndTime: 2019-10-10 10:20:52
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (878->38), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
	t80 = cos(pkin(6));
	t82 = sin(qJ(2));
	t86 = cos(qJ(1));
	t90 = t86 * t82;
	t83 = sin(qJ(1));
	t85 = cos(qJ(2));
	t91 = t83 * t85;
	t72 = t80 * t90 + t91;
	t78 = pkin(11) + qJ(4);
	t76 = sin(t78);
	t77 = cos(t78);
	t79 = sin(pkin(6));
	t93 = t79 * t86;
	t62 = t72 * t77 - t76 * t93;
	t96 = t79 * t82;
	t70 = t80 * t76 + t77 * t96;
	t60 = atan2(-t62, t70);
	t53 = sin(t60);
	t54 = cos(t60);
	t51 = -t53 * t62 + t54 * t70;
	t50 = 0.1e1 / t51 ^ 2;
	t89 = t86 * t85;
	t92 = t83 * t82;
	t74 = -t80 * t92 + t89;
	t95 = t79 * t83;
	t66 = t74 * t77 + t76 * t95;
	t103 = t50 * t66;
	t102 = t54 * t62;
	t65 = t74 * t76 - t77 * t95;
	t81 = sin(qJ(6));
	t73 = t80 * t91 + t90;
	t84 = cos(qJ(6));
	t97 = t73 * t84;
	t59 = t65 * t81 + t97;
	t56 = 0.1e1 / t59 ^ 2;
	t98 = t73 * t81;
	t58 = -t65 * t84 + t98;
	t101 = t56 * t58;
	t68 = 0.1e1 / t70 ^ 2;
	t100 = t62 * t68;
	t99 = t66 ^ 2 * t50;
	t94 = t79 * t85;
	t88 = t58 ^ 2 * t56 + 0.1e1;
	t87 = -t53 * t70 - t102;
	t61 = t72 * t76 + t77 * t93;
	t71 = t80 * t89 - t92;
	t69 = -t76 * t96 + t80 * t77;
	t67 = 0.1e1 / t70;
	t57 = 0.1e1 / (t62 ^ 2 * t68 + 0.1e1);
	t55 = 0.1e1 / t59;
	t52 = 0.1e1 / t88;
	t49 = 0.1e1 / t51;
	t48 = 0.1e1 / (0.1e1 + t99);
	t47 = (t94 * t100 - t67 * t71) * t77 * t57;
	t46 = (t69 * t100 + t61 * t67) * t57;
	t1 = [-t66 * t67 * t57, t47, 0, t46, 0, 0; (-t62 * t49 - (-t53 + (t67 * t102 + t53) * t57) * t99) * t48, (-t73 * t77 * t49 - ((-t53 * t71 + t54 * t94) * t77 + t87 * t47) * t103) * t48, 0, (-t65 * t49 - (t87 * t46 + t53 * t61 + t54 * t69) * t103) * t48, 0, 0; ((t61 * t84 + t71 * t81) * t55 - (-t61 * t81 + t71 * t84) * t101) * t52, ((t74 * t81 + t76 * t97) * t55 - (t74 * t84 - t76 * t98) * t101) * t52, 0, (-t81 * t101 - t84 * t55) * t66 * t52, 0, t88 * t52;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end