% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR7
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
%   Wie in S6RRPPRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPPRR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
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
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (58->16), mult. (119->35), div. (26->9), fcn. (180->9), ass. (0->26)
	t36 = cos(pkin(6));
	t35 = sin(pkin(6));
	t40 = cos(qJ(1));
	t43 = t40 * t35;
	t31 = atan2(-t43, -t36);
	t29 = sin(t31);
	t30 = cos(t31);
	t22 = -t29 * t43 - t30 * t36;
	t38 = sin(qJ(1));
	t48 = 0.1e1 / t22 ^ 2 * t38 ^ 2;
	t37 = sin(qJ(2));
	t42 = t40 * t37;
	t39 = cos(qJ(2));
	t44 = t38 * t39;
	t27 = t36 * t44 + t42;
	t25 = 0.1e1 / t27 ^ 2;
	t41 = t40 * t39;
	t45 = t38 * t37;
	t28 = -t36 * t45 + t41;
	t47 = t28 ^ 2 * t25;
	t33 = t35 ^ 2;
	t32 = 0.1e1 / (0.1e1 + t40 ^ 2 * t33 / t36 ^ 2);
	t46 = t32 / t36;
	t24 = 0.1e1 / t27;
	t23 = 0.1e1 / (0.1e1 + t47);
	t1 = [-t38 * t35 * t46, 0, 0, 0, 0, 0; (-0.1e1 / t22 * t43 + (t30 * t33 * t40 * t46 + (-t32 + 0.1e1) * t35 * t29) * t35 * t48) / (t33 * t48 + 0.1e1), 0, 0, 0, 0, 0; ((-t36 * t42 - t44) * t24 - (t36 * t41 - t45) * t28 * t25) * t23, (-t24 * t27 - t47) * t23, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (149->24), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->40)
	t52 = sin(qJ(2));
	t55 = cos(qJ(2));
	t56 = cos(qJ(1));
	t53 = sin(qJ(1));
	t60 = cos(pkin(6));
	t58 = t53 * t60;
	t43 = t56 * t52 + t55 * t58;
	t51 = sin(qJ(5));
	t54 = cos(qJ(5));
	t50 = sin(pkin(6));
	t62 = t50 * t53;
	t35 = t43 * t54 - t51 * t62;
	t33 = 0.1e1 / t35 ^ 2;
	t34 = t43 * t51 + t54 * t62;
	t67 = t33 * t34;
	t57 = t56 * t60;
	t41 = t52 * t57 + t53 * t55;
	t63 = t50 * t52;
	t39 = atan2(-t41, t63);
	t37 = cos(t39);
	t66 = t37 * t41;
	t36 = sin(t39);
	t30 = -t36 * t41 + t37 * t63;
	t29 = 0.1e1 / t30 ^ 2;
	t44 = -t52 * t58 + t56 * t55;
	t65 = t44 ^ 2 * t29;
	t47 = 0.1e1 / t50;
	t48 = 0.1e1 / t52;
	t64 = t47 * t48;
	t61 = t50 * t56;
	t59 = t34 ^ 2 * t33 + 0.1e1;
	t49 = 0.1e1 / t52 ^ 2;
	t40 = t53 * t52 - t55 * t57;
	t38 = 0.1e1 / (0.1e1 + t41 ^ 2 / t50 ^ 2 * t49);
	t32 = 0.1e1 / t35;
	t31 = 0.1e1 / t59;
	t28 = 0.1e1 / t30;
	t27 = 0.1e1 / (0.1e1 + t65);
	t26 = (t41 * t49 * t55 + t40 * t48) * t47 * t38;
	t1 = [-t44 * t38 * t64, t26, 0, 0, 0, 0; (-t41 * t28 - (-t36 + (t64 * t66 + t36) * t38) * t65) * t27, (-t43 * t28 - (t37 * t50 * t55 + t36 * t40 + (-t36 * t63 - t66) * t26) * t44 * t29) * t27, 0, 0, 0, 0; ((-t40 * t51 + t54 * t61) * t32 - (-t40 * t54 - t51 * t61) * t67) * t31, (t51 * t32 - t54 * t67) * t44 * t31, 0, 0, t59 * t31, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (459->37), mult. (1300->91), div. (85->9), fcn. (1858->13), ass. (0->54)
	t73 = cos(pkin(6));
	t80 = cos(qJ(2));
	t81 = cos(qJ(1));
	t84 = t81 * t80;
	t76 = sin(qJ(2));
	t77 = sin(qJ(1));
	t87 = t77 * t76;
	t66 = -t73 * t84 + t87;
	t75 = sin(qJ(5));
	t79 = cos(qJ(5));
	t72 = sin(pkin(6));
	t88 = t72 * t81;
	t56 = t66 * t75 - t79 * t88;
	t89 = t72 * t80;
	t64 = t73 * t79 - t75 * t89;
	t55 = atan2(-t56, t64);
	t52 = sin(t55);
	t53 = cos(t55);
	t46 = -t52 * t56 + t53 * t64;
	t45 = 0.1e1 / t46 ^ 2;
	t85 = t81 * t76;
	t86 = t77 * t80;
	t68 = t73 * t86 + t85;
	t90 = t72 * t77;
	t60 = t68 * t75 + t79 * t90;
	t97 = t45 * t60;
	t61 = t68 * t79 - t75 * t90;
	t69 = -t73 * t87 + t84;
	t74 = sin(qJ(6));
	t78 = cos(qJ(6));
	t51 = t61 * t78 + t69 * t74;
	t49 = 0.1e1 / t51 ^ 2;
	t50 = t61 * t74 - t69 * t78;
	t96 = t49 * t50;
	t95 = t53 * t56;
	t63 = 0.1e1 / t64 ^ 2;
	t94 = t56 * t63;
	t93 = t60 ^ 2 * t45;
	t92 = t69 * t79;
	t91 = t72 * t76;
	t83 = t50 ^ 2 * t49 + 0.1e1;
	t82 = -t52 * t64 - t95;
	t58 = t66 * t79 + t75 * t88;
	t67 = t73 * t85 + t86;
	t65 = -t73 * t75 - t79 * t89;
	t62 = 0.1e1 / t64;
	t54 = 0.1e1 / (t56 ^ 2 * t63 + 0.1e1);
	t48 = 0.1e1 / t51;
	t47 = 0.1e1 / t83;
	t44 = 0.1e1 / t46;
	t43 = 0.1e1 / (0.1e1 + t93);
	t42 = (-t62 * t67 + t91 * t94) * t75 * t54;
	t41 = (-t58 * t62 + t65 * t94) * t54;
	t1 = [-t60 * t62 * t54, t42, 0, 0, t41, 0; (-t56 * t44 - (-t52 + (t62 * t95 + t52) * t54) * t93) * t43, (t69 * t75 * t44 - ((-t52 * t67 + t53 * t91) * t75 + t82 * t42) * t97) * t43, 0, 0, (t61 * t44 - (t82 * t41 - t52 * t58 + t53 * t65) * t97) * t43, 0; ((-t58 * t74 + t67 * t78) * t48 - (-t58 * t78 - t67 * t74) * t96) * t47, ((t68 * t78 + t74 * t92) * t48 - (-t68 * t74 + t78 * t92) * t96) * t47, 0, 0, (-t74 * t48 + t78 * t96) * t60 * t47, t83 * t47;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end