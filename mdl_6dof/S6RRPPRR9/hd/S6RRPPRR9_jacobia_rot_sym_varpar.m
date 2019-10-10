% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR9
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
%   Wie in S6RRPPRR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR9_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR9_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
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
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
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
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (102->20), mult. (316->53), div. (70->11), fcn. (497->9), ass. (0->33)
	t41 = sin(qJ(2));
	t42 = sin(qJ(1));
	t43 = cos(qJ(2));
	t44 = cos(qJ(1));
	t47 = cos(pkin(6));
	t45 = t44 * t47;
	t29 = t41 * t45 + t42 * t43;
	t40 = sin(pkin(6));
	t48 = t40 * t41;
	t27 = atan2(-t29, t48);
	t24 = cos(t27);
	t52 = t24 * t29;
	t46 = t42 * t47;
	t31 = -t44 * t41 - t43 * t46;
	t34 = 0.1e1 / t40;
	t35 = 0.1e1 / t40 ^ 2;
	t39 = 0.1e1 / t42 ^ 2;
	t51 = 0.1e1 / (t31 ^ 2 * t39 * t35 + 0.1e1) * t34;
	t23 = sin(t27);
	t22 = -t23 * t29 + t24 * t48;
	t21 = 0.1e1 / t22 ^ 2;
	t32 = -t41 * t46 + t44 * t43;
	t50 = t32 ^ 2 * t21;
	t36 = 0.1e1 / t41;
	t49 = t34 * t36;
	t38 = 0.1e1 / t42;
	t37 = 0.1e1 / t41 ^ 2;
	t28 = t42 * t41 - t43 * t45;
	t25 = 0.1e1 / (t29 ^ 2 * t35 * t37 + 0.1e1);
	t20 = 0.1e1 / t22;
	t19 = 0.1e1 / (0.1e1 + t50);
	t18 = (t29 * t37 * t43 + t28 * t36) * t34 * t25;
	t1 = [-t32 * t25 * t49, t18, 0, 0, 0, 0; (-t29 * t20 - (-t23 + (t49 * t52 + t23) * t25) * t50) * t19, (t31 * t20 - (t24 * t40 * t43 + t23 * t28 + (-t23 * t48 - t52) * t18) * t32 * t21) * t19, 0, 0, 0, 0; (-t31 * t39 * t44 + t28 * t38) * t51, -t32 * t38 * t51, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (149->22), mult. (451->62), div. (72->11), fcn. (673->11), ass. (0->42)
	t52 = cos(pkin(6));
	t57 = cos(qJ(2));
	t58 = cos(qJ(1));
	t60 = t58 * t57;
	t54 = sin(qJ(2));
	t55 = sin(qJ(1));
	t63 = t55 * t54;
	t46 = -t52 * t63 + t60;
	t53 = sin(qJ(5));
	t56 = cos(qJ(5));
	t51 = sin(pkin(6));
	t66 = t51 * t55;
	t37 = t46 * t53 + t56 * t66;
	t35 = 0.1e1 / t37 ^ 2;
	t36 = -t46 * t56 + t53 * t66;
	t70 = t35 * t36;
	t42 = -t52 * t60 + t63;
	t65 = t51 * t57;
	t41 = atan2(t42, t65);
	t39 = cos(t41);
	t69 = t39 * t42;
	t38 = sin(t41);
	t32 = t38 * t42 + t39 * t65;
	t31 = 0.1e1 / t32 ^ 2;
	t61 = t58 * t54;
	t62 = t55 * t57;
	t44 = t52 * t62 + t61;
	t68 = t44 ^ 2 * t31;
	t48 = 0.1e1 / t51;
	t49 = 0.1e1 / t57;
	t67 = t48 * t49;
	t64 = t51 * t58;
	t59 = t36 ^ 2 * t35 + 0.1e1;
	t50 = 0.1e1 / t57 ^ 2;
	t43 = t52 * t61 + t62;
	t40 = 0.1e1 / (0.1e1 + t42 ^ 2 / t51 ^ 2 * t50);
	t34 = 0.1e1 / t37;
	t33 = 0.1e1 / t59;
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / (0.1e1 + t68);
	t28 = (t42 * t50 * t54 + t43 * t49) * t48 * t40;
	t1 = [t44 * t40 * t67, t28, 0, 0, 0, 0; (t42 * t30 + (t38 + (t67 * t69 - t38) * t40) * t68) * t29, (-t46 * t30 + (-t39 * t51 * t54 + t38 * t43 + (-t38 * t65 + t69) * t28) * t44 * t31) * t29, 0, 0, 0, 0; ((t43 * t56 + t53 * t64) * t34 - (-t43 * t53 + t56 * t64) * t70) * t33, (t56 * t34 + t53 * t70) * t44 * t33, 0, 0, t59 * t33, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (459->36), mult. (1300->90), div. (85->9), fcn. (1858->13), ass. (0->55)
	t75 = cos(pkin(6));
	t78 = sin(qJ(2));
	t83 = cos(qJ(1));
	t89 = t83 * t78;
	t79 = sin(qJ(1));
	t82 = cos(qJ(2));
	t90 = t79 * t82;
	t70 = t75 * t89 + t90;
	t77 = sin(qJ(5));
	t81 = cos(qJ(5));
	t74 = sin(pkin(6));
	t92 = t74 * t83;
	t61 = t70 * t81 + t77 * t92;
	t94 = t74 * t81;
	t67 = t75 * t77 - t78 * t94;
	t58 = atan2(t61, t67);
	t55 = sin(t58);
	t56 = cos(t58);
	t49 = t55 * t61 + t56 * t67;
	t48 = 0.1e1 / t49 ^ 2;
	t88 = t83 * t82;
	t91 = t78 * t79;
	t86 = -t75 * t91 + t88;
	t95 = t74 * t77;
	t59 = t79 * t95 - t86 * t81;
	t102 = t48 * t59;
	t101 = t48 * t59 ^ 2;
	t60 = t86 * t77 + t79 * t94;
	t80 = cos(qJ(6));
	t71 = -t75 * t90 - t89;
	t76 = sin(qJ(6));
	t97 = t71 * t76;
	t54 = t60 * t80 + t97;
	t52 = 0.1e1 / t54 ^ 2;
	t96 = t71 * t80;
	t53 = t60 * t76 - t96;
	t100 = t52 * t53;
	t99 = t56 * t61;
	t66 = 0.1e1 / t67 ^ 2;
	t98 = t61 * t66;
	t93 = t74 * t82;
	t87 = t52 * t53 ^ 2 + 0.1e1;
	t85 = -t55 * t67 + t99;
	t84 = -t70 * t77 + t81 * t92;
	t69 = -t75 * t88 + t91;
	t68 = t75 * t81 + t78 * t95;
	t65 = 0.1e1 / t67;
	t57 = 0.1e1 / (t61 ^ 2 * t66 + 0.1e1);
	t51 = 0.1e1 / t54;
	t50 = 0.1e1 / t87;
	t47 = 0.1e1 / t49;
	t46 = 0.1e1 / (0.1e1 + t101);
	t45 = (-t65 * t69 + t93 * t98) * t81 * t57;
	t44 = (t65 * t84 - t68 * t98) * t57;
	t1 = [-t59 * t65 * t57, t45, 0, 0, t44, 0; (t61 * t47 - (-t55 + (-t65 * t99 + t55) * t57) * t101) * t46, (-t71 * t81 * t47 - ((-t55 * t69 - t56 * t93) * t81 + t85 * t45) * t102) * t46, 0, 0, (t60 * t47 - (t85 * t44 + t55 * t84 + t56 * t68) * t102) * t46, 0; ((-t69 * t80 + t76 * t84) * t51 - (t69 * t76 + t80 * t84) * t100) * t50, ((t77 * t97 + t86 * t80) * t51 - (-t86 * t76 + t77 * t96) * t100) * t50, 0, 0, (t80 * t100 - t76 * t51) * t59 * t50, t87 * t50;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end