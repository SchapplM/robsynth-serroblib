% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR5
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
%   Wie in S6RRPPRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:43
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
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
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (128->20), mult. (330->51), div. (64->11), fcn. (506->9), ass. (0->33)
	t40 = sin(qJ(2));
	t41 = sin(qJ(1));
	t42 = cos(qJ(2));
	t43 = cos(qJ(1));
	t47 = cos(pkin(6));
	t45 = t43 * t47;
	t30 = t41 * t40 - t42 * t45;
	t39 = sin(pkin(6));
	t48 = t39 * t42;
	t26 = atan2(-t30, -t48);
	t25 = cos(t26);
	t53 = t25 * t30;
	t46 = t41 * t47;
	t34 = -t40 * t46 + t43 * t42;
	t29 = 0.1e1 / t34 ^ 2;
	t44 = t39 ^ 2;
	t52 = 0.1e1 / (t41 ^ 2 * t44 * t29 + 0.1e1) * t39;
	t51 = t29 * t41;
	t24 = sin(t26);
	t23 = -t24 * t30 - t25 * t48;
	t22 = 0.1e1 / t23 ^ 2;
	t33 = t43 * t40 + t42 * t46;
	t50 = t33 ^ 2 * t22;
	t36 = 0.1e1 / t39;
	t37 = 0.1e1 / t42;
	t49 = t36 * t37;
	t38 = 0.1e1 / t42 ^ 2;
	t32 = t40 * t45 + t41 * t42;
	t28 = 0.1e1 / (0.1e1 + t30 ^ 2 / t44 * t38);
	t21 = 0.1e1 / t23;
	t20 = 0.1e1 / (0.1e1 + t50);
	t19 = (t30 * t38 * t40 + t32 * t37) * t36 * t28;
	t1 = [t33 * t28 * t49, t19, 0, 0, 0, 0; (-t30 * t21 - (-t24 + (-t49 * t53 + t24) * t28) * t50) * t20, (t34 * t21 - (t25 * t39 * t40 - t24 * t32 + (t24 * t48 - t53) * t19) * t33 * t22) * t20, 0, 0, 0, 0; (-t43 / t34 - t32 * t51) * t52, -t33 * t51 * t52, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (58->16), mult. (119->35), div. (26->9), fcn. (180->9), ass. (0->26)
	t35 = cos(pkin(6));
	t34 = sin(pkin(6));
	t39 = cos(qJ(1));
	t42 = t39 * t34;
	t30 = atan2(-t42, -t35);
	t28 = sin(t30);
	t29 = cos(t30);
	t21 = -t28 * t42 - t29 * t35;
	t37 = sin(qJ(1));
	t47 = 0.1e1 / t21 ^ 2 * t37 ^ 2;
	t38 = cos(qJ(2));
	t40 = t39 * t38;
	t36 = sin(qJ(2));
	t44 = t37 * t36;
	t27 = -t35 * t44 + t40;
	t25 = 0.1e1 / t27 ^ 2;
	t41 = t39 * t36;
	t43 = t37 * t38;
	t26 = -t35 * t43 - t41;
	t46 = t26 ^ 2 * t25;
	t32 = t34 ^ 2;
	t31 = 0.1e1 / (0.1e1 + t39 ^ 2 * t32 / t35 ^ 2);
	t45 = t31 / t35;
	t24 = 0.1e1 / t27;
	t22 = 0.1e1 / (0.1e1 + t46);
	t1 = [-t37 * t34 * t45, 0, 0, 0, 0, 0; (-0.1e1 / t21 * t42 + (t29 * t32 * t39 * t45 + (-t31 + 0.1e1) * t34 * t28) * t34 * t47) / (t32 * t47 + 0.1e1), 0, 0, 0, 0, 0; ((-t35 * t40 + t44) * t24 - (-t35 * t41 - t43) * t26 * t25) * t22, (-t24 * t27 - t46) * t22, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (149->22), mult. (451->62), div. (72->11), fcn. (673->11), ass. (0->42)
	t50 = cos(pkin(6));
	t55 = cos(qJ(2));
	t56 = cos(qJ(1));
	t58 = t56 * t55;
	t52 = sin(qJ(2));
	t53 = sin(qJ(1));
	t61 = t53 * t52;
	t44 = -t50 * t61 + t58;
	t51 = sin(qJ(5));
	t54 = cos(qJ(5));
	t49 = sin(pkin(6));
	t64 = t49 * t53;
	t35 = t44 * t54 - t51 * t64;
	t33 = 0.1e1 / t35 ^ 2;
	t34 = t44 * t51 + t54 * t64;
	t68 = t33 * t34;
	t40 = -t50 * t58 + t61;
	t63 = t49 * t55;
	t39 = atan2(t40, t63);
	t37 = cos(t39);
	t67 = t37 * t40;
	t36 = sin(t39);
	t30 = t36 * t40 + t37 * t63;
	t29 = 0.1e1 / t30 ^ 2;
	t59 = t56 * t52;
	t60 = t53 * t55;
	t42 = t50 * t60 + t59;
	t66 = t42 ^ 2 * t29;
	t46 = 0.1e1 / t49;
	t47 = 0.1e1 / t55;
	t65 = t46 * t47;
	t62 = t49 * t56;
	t57 = t34 ^ 2 * t33 + 0.1e1;
	t48 = 0.1e1 / t55 ^ 2;
	t41 = t50 * t59 + t60;
	t38 = 0.1e1 / (0.1e1 + t40 ^ 2 / t49 ^ 2 * t48);
	t32 = 0.1e1 / t35;
	t31 = 0.1e1 / t57;
	t28 = 0.1e1 / t30;
	t27 = 0.1e1 / (0.1e1 + t66);
	t26 = (t40 * t48 * t52 + t41 * t47) * t46 * t38;
	t1 = [t42 * t38 * t65, t26, 0, 0, 0, 0; (t40 * t28 + (t36 + (t65 * t67 - t36) * t38) * t66) * t27, (-t44 * t28 + (-t37 * t49 * t52 + t36 * t41 + (-t36 * t63 + t67) * t26) * t42 * t29) * t27, 0, 0, 0, 0; ((-t41 * t51 + t54 * t62) * t32 - (-t41 * t54 - t51 * t62) * t68) * t31, (-t51 * t32 + t54 * t68) * t42 * t31, 0, 0, t57 * t31, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (459->37), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->54)
	t75 = cos(pkin(6));
	t78 = sin(qJ(2));
	t83 = cos(qJ(1));
	t87 = t83 * t78;
	t79 = sin(qJ(1));
	t82 = cos(qJ(2));
	t88 = t79 * t82;
	t68 = t75 * t87 + t88;
	t77 = sin(qJ(5));
	t81 = cos(qJ(5));
	t74 = sin(pkin(6));
	t90 = t74 * t83;
	t57 = t68 * t77 - t81 * t90;
	t93 = t74 * t77;
	t65 = t75 * t81 + t78 * t93;
	t56 = atan2(-t57, t65);
	t53 = sin(t56);
	t54 = cos(t56);
	t47 = -t53 * t57 + t54 * t65;
	t46 = 0.1e1 / t47 ^ 2;
	t86 = t83 * t82;
	t89 = t79 * t78;
	t70 = -t75 * t89 + t86;
	t92 = t74 * t81;
	t61 = t70 * t77 + t79 * t92;
	t99 = t46 * t61;
	t62 = t70 * t81 - t79 * t93;
	t69 = -t75 * t88 - t87;
	t76 = sin(qJ(6));
	t80 = cos(qJ(6));
	t52 = t62 * t80 + t69 * t76;
	t50 = 0.1e1 / t52 ^ 2;
	t51 = t62 * t76 - t69 * t80;
	t98 = t50 * t51;
	t97 = t54 * t57;
	t64 = 0.1e1 / t65 ^ 2;
	t96 = t57 * t64;
	t95 = t61 ^ 2 * t46;
	t94 = t69 * t81;
	t91 = t74 * t82;
	t85 = t51 ^ 2 * t50 + 0.1e1;
	t84 = -t53 * t65 - t97;
	t59 = t68 * t81 + t77 * t90;
	t67 = -t75 * t86 + t89;
	t66 = -t75 * t77 + t78 * t92;
	t63 = 0.1e1 / t65;
	t55 = 0.1e1 / (t57 ^ 2 * t64 + 0.1e1);
	t49 = 0.1e1 / t52;
	t48 = 0.1e1 / t85;
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / (0.1e1 + t95);
	t43 = (t63 * t67 + t91 * t96) * t77 * t55;
	t42 = (-t59 * t63 + t66 * t96) * t55;
	t1 = [-t61 * t63 * t55, t43, 0, 0, t42, 0; (-t57 * t45 - (-t53 + (t63 * t97 + t53) * t55) * t95) * t44, (t69 * t77 * t45 - ((t53 * t67 + t54 * t91) * t77 + t84 * t43) * t99) * t44, 0, 0, (t62 * t45 - (t42 * t84 - t53 * t59 + t54 * t66) * t99) * t44, 0; ((-t59 * t76 - t67 * t80) * t49 - (-t59 * t80 + t67 * t76) * t98) * t48, ((t70 * t80 + t76 * t94) * t49 - (-t70 * t76 + t80 * t94) * t98) * t48, 0, 0, (-t76 * t49 + t80 * t98) * t61 * t48, t85 * t48;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end