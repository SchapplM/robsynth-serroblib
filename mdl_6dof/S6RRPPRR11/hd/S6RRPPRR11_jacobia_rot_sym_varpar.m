% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR11
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
%   Wie in S6RRPPRR11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:53
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR11_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR11_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
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
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (126->20), mult. (316->53), div. (70->11), fcn. (497->9), ass. (0->33)
	t43 = sin(qJ(2));
	t44 = sin(qJ(1));
	t45 = cos(qJ(2));
	t46 = cos(qJ(1));
	t49 = cos(pkin(6));
	t47 = t46 * t49;
	t30 = t44 * t43 - t45 * t47;
	t42 = sin(pkin(6));
	t50 = t42 * t45;
	t27 = atan2(-t30, -t50);
	t26 = cos(t27);
	t54 = t26 * t30;
	t48 = t44 * t49;
	t34 = -t43 * t48 + t46 * t45;
	t36 = 0.1e1 / t42;
	t37 = 0.1e1 / t42 ^ 2;
	t39 = 0.1e1 / t44 ^ 2;
	t53 = 0.1e1 / (t34 ^ 2 * t39 * t37 + 0.1e1) * t36;
	t25 = sin(t27);
	t24 = -t25 * t30 - t26 * t50;
	t23 = 0.1e1 / t24 ^ 2;
	t33 = t46 * t43 + t45 * t48;
	t52 = t33 ^ 2 * t23;
	t40 = 0.1e1 / t45;
	t51 = t36 * t40;
	t41 = 0.1e1 / t45 ^ 2;
	t38 = 0.1e1 / t44;
	t32 = t43 * t47 + t44 * t45;
	t28 = 0.1e1 / (t30 ^ 2 * t37 * t41 + 0.1e1);
	t22 = 0.1e1 / t24;
	t21 = 0.1e1 / (0.1e1 + t52);
	t20 = (t30 * t41 * t43 + t32 * t40) * t36 * t28;
	t1 = [t33 * t28 * t51, t20, 0, 0, 0, 0; (-t30 * t22 - (-t25 + (-t51 * t54 + t25) * t28) * t52) * t21, (t34 * t22 - (t26 * t42 * t43 - t25 * t32 + (t25 * t50 - t54) * t20) * t33 * t23) * t21, 0, 0, 0, 0; (-t34 * t39 * t46 - t32 * t38) * t53, -t33 * t38 * t53, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (132->24), mult. (403->63), div. (67->11), fcn. (610->11), ass. (0->39)
	t47 = sin(qJ(2));
	t49 = cos(qJ(2));
	t50 = cos(qJ(1));
	t48 = sin(qJ(1));
	t53 = cos(pkin(6));
	t52 = t48 * t53;
	t37 = t50 * t47 + t49 * t52;
	t44 = sin(pkin(11));
	t46 = cos(pkin(11));
	t45 = sin(pkin(6));
	t55 = t45 * t48;
	t29 = t37 * t44 + t46 * t55;
	t27 = 0.1e1 / t29 ^ 2;
	t28 = -t37 * t46 + t44 * t55;
	t60 = t27 * t28;
	t51 = t50 * t53;
	t35 = t47 * t51 + t48 * t49;
	t56 = t45 * t47;
	t33 = atan2(-t35, t56);
	t31 = cos(t33);
	t59 = t31 * t35;
	t30 = sin(t33);
	t24 = -t30 * t35 + t31 * t56;
	t23 = 0.1e1 / t24 ^ 2;
	t38 = -t47 * t52 + t50 * t49;
	t58 = t38 ^ 2 * t23;
	t41 = 0.1e1 / t45;
	t42 = 0.1e1 / t47;
	t57 = t41 * t42;
	t54 = t45 * t50;
	t43 = 0.1e1 / t47 ^ 2;
	t34 = t48 * t47 - t49 * t51;
	t32 = 0.1e1 / (0.1e1 + t35 ^ 2 / t45 ^ 2 * t43);
	t26 = 0.1e1 / t29;
	t25 = 0.1e1 / (t28 ^ 2 * t27 + 0.1e1);
	t22 = 0.1e1 / t24;
	t21 = 0.1e1 / (0.1e1 + t58);
	t20 = (t35 * t43 * t49 + t34 * t42) * t41 * t32;
	t1 = [-t38 * t32 * t57, t20, 0, 0, 0, 0; (-t35 * t22 - (-t30 + (t57 * t59 + t30) * t32) * t58) * t21, (-t37 * t22 - (t31 * t45 * t49 + t30 * t34 + (-t30 * t56 - t59) * t20) * t38 * t23) * t21, 0, 0, 0, 0; ((t34 * t46 + t44 * t54) * t26 - (-t34 * t44 + t46 * t54) * t60) * t25, (-t46 * t26 - t44 * t60) * t38 * t25, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (197->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
	t59 = sin(qJ(2));
	t61 = cos(qJ(2));
	t62 = cos(qJ(1));
	t60 = sin(qJ(1));
	t66 = cos(pkin(6));
	t64 = t60 * t66;
	t48 = t62 * t59 + t61 * t64;
	t55 = pkin(11) + qJ(5);
	t52 = sin(t55);
	t53 = cos(t55);
	t58 = sin(pkin(6));
	t68 = t58 * t60;
	t40 = t48 * t52 + t53 * t68;
	t38 = 0.1e1 / t40 ^ 2;
	t39 = -t48 * t53 + t52 * t68;
	t73 = t38 * t39;
	t63 = t62 * t66;
	t46 = t59 * t63 + t60 * t61;
	t69 = t58 * t59;
	t44 = atan2(-t46, t69);
	t42 = cos(t44);
	t72 = t42 * t46;
	t41 = sin(t44);
	t35 = -t41 * t46 + t42 * t69;
	t34 = 0.1e1 / t35 ^ 2;
	t49 = -t59 * t64 + t62 * t61;
	t71 = t49 ^ 2 * t34;
	t54 = 0.1e1 / t58;
	t56 = 0.1e1 / t59;
	t70 = t54 * t56;
	t67 = t58 * t62;
	t65 = t39 ^ 2 * t38 + 0.1e1;
	t57 = 0.1e1 / t59 ^ 2;
	t45 = t60 * t59 - t61 * t63;
	t43 = 0.1e1 / (0.1e1 + t46 ^ 2 / t58 ^ 2 * t57);
	t37 = 0.1e1 / t40;
	t36 = 0.1e1 / t65;
	t33 = 0.1e1 / t35;
	t32 = 0.1e1 / (0.1e1 + t71);
	t31 = (t46 * t57 * t61 + t45 * t56) * t54 * t43;
	t1 = [-t49 * t43 * t70, t31, 0, 0, 0, 0; (-t46 * t33 - (-t41 + (t70 * t72 + t41) * t43) * t71) * t32, (-t48 * t33 - (t42 * t58 * t61 + t41 * t45 + (-t41 * t69 - t72) * t31) * t49 * t34) * t32, 0, 0, 0, 0; ((t45 * t53 + t52 * t67) * t37 - (-t45 * t52 + t53 * t67) * t73) * t36, (-t53 * t37 - t52 * t73) * t49 * t36, 0, 0, t65 * t36, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (878->37), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->56)
	t80 = cos(pkin(6));
	t85 = cos(qJ(2));
	t86 = cos(qJ(1));
	t91 = t86 * t85;
	t82 = sin(qJ(2));
	t83 = sin(qJ(1));
	t94 = t83 * t82;
	t72 = -t80 * t91 + t94;
	t78 = pkin(11) + qJ(5);
	t76 = sin(t78);
	t77 = cos(t78);
	t79 = sin(pkin(6));
	t95 = t79 * t86;
	t64 = t72 * t77 + t76 * t95;
	t96 = t79 * t85;
	t70 = t80 * t76 + t77 * t96;
	t61 = atan2(t64, t70);
	t54 = sin(t61);
	t55 = cos(t61);
	t52 = t54 * t64 + t55 * t70;
	t51 = 0.1e1 / t52 ^ 2;
	t92 = t86 * t82;
	t93 = t83 * t85;
	t87 = t80 * t93 + t92;
	t97 = t79 * t83;
	t62 = t76 * t97 - t87 * t77;
	t105 = t51 * t62;
	t104 = t55 * t64;
	t74 = -t80 * t94 + t91;
	t81 = sin(qJ(6));
	t100 = t74 * t81;
	t63 = t87 * t76 + t77 * t97;
	t84 = cos(qJ(6));
	t60 = t63 * t84 + t100;
	t57 = 0.1e1 / t60 ^ 2;
	t99 = t74 * t84;
	t59 = t63 * t81 - t99;
	t103 = t57 * t59;
	t102 = t62 ^ 2 * t51;
	t69 = 0.1e1 / t70 ^ 2;
	t101 = t64 * t69;
	t98 = t79 * t82;
	t90 = t59 ^ 2 * t57 + 0.1e1;
	t89 = -t54 * t70 + t104;
	t88 = -t72 * t76 + t77 * t95;
	t73 = t80 * t92 + t93;
	t71 = -t76 * t96 + t80 * t77;
	t68 = 0.1e1 / t70;
	t58 = 0.1e1 / (t64 ^ 2 * t69 + 0.1e1);
	t56 = 0.1e1 / t60;
	t53 = 0.1e1 / t90;
	t50 = 0.1e1 / t52;
	t49 = 0.1e1 / (0.1e1 + t102);
	t48 = (t98 * t101 + t68 * t73) * t77 * t58;
	t47 = (-t71 * t101 + t68 * t88) * t58;
	t1 = [-t62 * t68 * t58, t48, 0, 0, t47, 0; (t64 * t50 - (-t54 + (-t68 * t104 + t54) * t58) * t102) * t49, (-t74 * t77 * t50 - ((t54 * t73 - t55 * t98) * t77 + t89 * t48) * t105) * t49, 0, 0, (t63 * t50 - (t89 * t47 + t54 * t88 + t55 * t71) * t105) * t49, 0; ((t73 * t84 + t81 * t88) * t56 - (-t73 * t81 + t84 * t88) * t103) * t53, ((t76 * t100 + t87 * t84) * t56 - (t76 * t99 - t87 * t81) * t103) * t53, 0, 0, (t84 * t103 - t81 * t56) * t62 * t53, t90 * t53;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end