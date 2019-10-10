% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR2
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
%   Wie in S6PRRPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:09
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
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
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (131->19), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->30)
	t50 = sin(pkin(10));
	t51 = sin(pkin(6));
	t60 = t50 * t51;
	t55 = cos(qJ(2));
	t59 = t51 * t55;
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t58 = t53 * t54;
	t57 = t53 * t55;
	t52 = cos(pkin(10));
	t43 = -t50 * t58 + t52 * t55;
	t48 = qJ(3) + pkin(11);
	t45 = sin(t48);
	t46 = cos(t48);
	t34 = t43 * t46 + t45 * t60;
	t32 = 0.1e1 / t34 ^ 2;
	t33 = t43 * t45 - t46 * t60;
	t56 = t33 ^ 2 * t32 + 0.1e1;
	t49 = 0.1e1 / t55 ^ 2;
	t42 = t50 * t57 + t52 * t54;
	t41 = t50 * t55 + t52 * t58;
	t39 = t50 * t54 - t52 * t57;
	t37 = atan2(-t39, -t59);
	t36 = cos(t37);
	t35 = sin(t37);
	t31 = 0.1e1 / t56;
	t30 = -t35 * t39 - t36 * t59;
	t29 = 0.1e1 / t30 ^ 2;
	t27 = (t41 / t55 + t54 * t39 * t49) / t51 / (0.1e1 + t39 ^ 2 / t51 ^ 2 * t49);
	t1 = [0, t27, 0, 0, 0, 0; 0, (t43 / t30 - (t36 * t51 * t54 - t35 * t41 + (t35 * t59 - t36 * t39) * t27) * t42 * t29) / (t42 ^ 2 * t29 + 0.1e1), 0, 0, 0, 0; 0, (-t45 / t34 + t46 * t33 * t32) * t42 * t31, t56 * t31, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (523->27), mult. (731->69), div. (57->9), fcn. (1053->11), ass. (0->44)
	t65 = sin(pkin(10));
	t67 = cos(pkin(10));
	t70 = cos(qJ(2));
	t68 = cos(pkin(6));
	t69 = sin(qJ(2));
	t73 = t68 * t69;
	t59 = t65 * t70 + t67 * t73;
	t64 = qJ(3) + pkin(11);
	t62 = sin(t64);
	t63 = cos(t64);
	t66 = sin(pkin(6));
	t76 = t66 * t67;
	t47 = t59 * t62 + t63 * t76;
	t75 = t66 * t69;
	t54 = t62 * t75 - t68 * t63;
	t45 = atan2(-t47, t54);
	t42 = sin(t45);
	t43 = cos(t45);
	t41 = -t42 * t47 + t43 * t54;
	t40 = 0.1e1 / t41 ^ 2;
	t61 = -t65 * t73 + t67 * t70;
	t77 = t65 * t66;
	t50 = t61 * t62 - t63 * t77;
	t79 = t40 * t50;
	t53 = 0.1e1 / t54 ^ 2;
	t78 = t47 * t53;
	t74 = t66 * t70;
	t72 = t68 * t70;
	t71 = -t42 * t54 - t43 * t47;
	t60 = t65 * t72 + t67 * t69;
	t58 = -t65 * t69 + t67 * t72;
	t57 = 0.1e1 / t60 ^ 2;
	t56 = 0.1e1 / t60;
	t55 = t68 * t62 + t63 * t75;
	t52 = 0.1e1 / t54;
	t51 = t61 * t63 + t62 * t77;
	t49 = t59 * t63 - t62 * t76;
	t46 = 0.1e1 / (t51 ^ 2 * t57 + 0.1e1);
	t44 = 0.1e1 / (t47 ^ 2 * t53 + 0.1e1);
	t39 = 0.1e1 / t41;
	t38 = 0.1e1 / (t50 ^ 2 * t40 + 0.1e1);
	t37 = (-t52 * t58 + t74 * t78) * t62 * t44;
	t36 = (-t49 * t52 + t55 * t78) * t44;
	t1 = [0, t37, t36, 0, 0, 0; 0, (-t60 * t62 * t39 - ((-t42 * t58 + t43 * t74) * t62 + t71 * t37) * t79) * t38, (t51 * t39 - (t71 * t36 - t42 * t49 + t43 * t55) * t79) * t38, 0, 0, 0; 0, (-t51 * t57 * t61 - t56 * t60 * t63) * t46, -t50 * t56 * t46, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (641->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->52)
	t72 = sin(pkin(10));
	t74 = cos(pkin(10));
	t79 = cos(qJ(2));
	t75 = cos(pkin(6));
	t77 = sin(qJ(2));
	t83 = t75 * t77;
	t65 = t72 * t79 + t74 * t83;
	t71 = qJ(3) + pkin(11);
	t69 = sin(t71);
	t70 = cos(t71);
	t73 = sin(pkin(6));
	t86 = t73 * t74;
	t56 = t65 * t70 - t69 * t86;
	t85 = t73 * t77;
	t63 = t69 * t75 + t70 * t85;
	t54 = atan2(-t56, t63);
	t49 = sin(t54);
	t50 = cos(t54);
	t45 = -t49 * t56 + t50 * t63;
	t44 = 0.1e1 / t45 ^ 2;
	t67 = -t72 * t83 + t74 * t79;
	t87 = t72 * t73;
	t59 = t67 * t70 + t69 * t87;
	t92 = t44 * t59;
	t58 = t67 * t69 - t70 * t87;
	t76 = sin(qJ(6));
	t82 = t75 * t79;
	t66 = t72 * t82 + t74 * t77;
	t78 = cos(qJ(6));
	t88 = t66 * t78;
	t52 = t58 * t76 + t88;
	t48 = 0.1e1 / t52 ^ 2;
	t89 = t66 * t76;
	t51 = -t58 * t78 + t89;
	t91 = t48 * t51;
	t61 = 0.1e1 / t63 ^ 2;
	t90 = t56 * t61;
	t84 = t73 * t79;
	t81 = t48 * t51 ^ 2 + 0.1e1;
	t80 = -t49 * t63 - t50 * t56;
	t64 = -t72 * t77 + t74 * t82;
	t62 = -t69 * t85 + t70 * t75;
	t60 = 0.1e1 / t63;
	t55 = t65 * t69 + t70 * t86;
	t53 = 0.1e1 / (t56 ^ 2 * t61 + 0.1e1);
	t47 = 0.1e1 / t52;
	t46 = 0.1e1 / t81;
	t43 = 0.1e1 / t45;
	t42 = 0.1e1 / (t44 * t59 ^ 2 + 0.1e1);
	t41 = (-t60 * t64 + t84 * t90) * t70 * t53;
	t40 = (t55 * t60 + t62 * t90) * t53;
	t1 = [0, t41, t40, 0, 0, 0; 0, (-t66 * t70 * t43 - ((-t49 * t64 + t50 * t84) * t70 + t80 * t41) * t92) * t42, (-t58 * t43 - (t80 * t40 + t49 * t55 + t50 * t62) * t92) * t42, 0, 0, 0; 0, ((t67 * t76 + t69 * t88) * t47 - (t67 * t78 - t69 * t89) * t91) * t46, (-t47 * t78 - t76 * t91) * t59 * t46, 0, 0, t81 * t46;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end