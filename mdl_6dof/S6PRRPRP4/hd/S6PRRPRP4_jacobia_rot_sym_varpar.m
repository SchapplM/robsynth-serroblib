% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRP4
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
%   Wie in S6PRRPRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:21
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRP4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
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
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (253->26), mult. (731->68), div. (57->9), fcn. (1053->11), ass. (0->42)
	t53 = sin(pkin(10));
	t55 = cos(pkin(10));
	t60 = cos(qJ(2));
	t56 = cos(pkin(6));
	t58 = sin(qJ(2));
	t63 = t56 * t58;
	t48 = t53 * t60 + t55 * t63;
	t57 = sin(qJ(3));
	t54 = sin(pkin(6));
	t59 = cos(qJ(3));
	t65 = t54 * t59;
	t38 = t48 * t57 + t55 * t65;
	t66 = t54 * t57;
	t51 = -t56 * t59 + t58 * t66;
	t37 = atan2(-t38, t51);
	t33 = sin(t37);
	t34 = cos(t37);
	t32 = -t33 * t38 + t34 * t51;
	t31 = 0.1e1 / t32 ^ 2;
	t50 = -t53 * t63 + t55 * t60;
	t41 = t50 * t57 - t53 * t65;
	t68 = t31 * t41;
	t46 = 0.1e1 / t51 ^ 2;
	t67 = t38 * t46;
	t64 = t54 * t60;
	t62 = t56 * t60;
	t61 = -t33 * t51 - t34 * t38;
	t52 = t56 * t57 + t58 * t65;
	t49 = t53 * t62 + t55 * t58;
	t47 = -t53 * t58 + t55 * t62;
	t45 = 0.1e1 / t51;
	t44 = 0.1e1 / t49 ^ 2;
	t43 = 0.1e1 / t49;
	t42 = t50 * t59 + t53 * t66;
	t40 = t48 * t59 - t55 * t66;
	t36 = 0.1e1 / (t38 ^ 2 * t46 + 0.1e1);
	t35 = 0.1e1 / (t42 ^ 2 * t44 + 0.1e1);
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / (t41 ^ 2 * t31 + 0.1e1);
	t28 = (-t45 * t47 + t64 * t67) * t57 * t36;
	t27 = (-t40 * t45 + t52 * t67) * t36;
	t1 = [0, t28, t27, 0, 0, 0; 0, (-t49 * t57 * t30 - ((-t33 * t47 + t34 * t64) * t57 + t61 * t28) * t68) * t29, (t42 * t30 - (t61 * t27 - t33 * t40 + t34 * t52) * t68) * t29, 0, 0, 0; 0, (-t42 * t44 * t50 - t43 * t49 * t59) * t35, -t41 * t43 * t35, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:55
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (334->30), mult. (949->78), div. (65->9), fcn. (1349->13), ass. (0->50)
	t64 = sin(pkin(10));
	t66 = cos(pkin(10));
	t73 = cos(qJ(2));
	t67 = cos(pkin(6));
	t70 = sin(qJ(2));
	t77 = t67 * t70;
	t58 = t64 * t73 + t66 * t77;
	t72 = cos(qJ(3));
	t65 = sin(pkin(6));
	t69 = sin(qJ(3));
	t80 = t65 * t69;
	t51 = t58 * t72 - t66 * t80;
	t79 = t65 * t72;
	t62 = t67 * t69 + t70 * t79;
	t49 = atan2(-t51, t62);
	t46 = sin(t49);
	t47 = cos(t49);
	t40 = -t46 * t51 + t47 * t62;
	t39 = 0.1e1 / t40 ^ 2;
	t60 = -t64 * t77 + t66 * t73;
	t54 = t60 * t72 + t64 * t80;
	t85 = t39 * t54;
	t53 = t60 * t69 - t64 * t79;
	t68 = sin(qJ(5));
	t76 = t67 * t73;
	t59 = t64 * t76 + t66 * t70;
	t71 = cos(qJ(5));
	t81 = t59 * t71;
	t45 = t53 * t68 + t81;
	t43 = 0.1e1 / t45 ^ 2;
	t82 = t59 * t68;
	t44 = -t53 * t71 + t82;
	t84 = t43 * t44;
	t56 = 0.1e1 / t62 ^ 2;
	t83 = t51 * t56;
	t78 = t65 * t73;
	t75 = t44 ^ 2 * t43 + 0.1e1;
	t74 = -t46 * t62 - t47 * t51;
	t61 = t67 * t72 - t70 * t80;
	t57 = -t64 * t70 + t66 * t76;
	t55 = 0.1e1 / t62;
	t50 = t58 * t69 + t66 * t79;
	t48 = 0.1e1 / (t51 ^ 2 * t56 + 0.1e1);
	t42 = 0.1e1 / t45;
	t41 = 0.1e1 / t75;
	t38 = 0.1e1 / t40;
	t37 = 0.1e1 / (t54 ^ 2 * t39 + 0.1e1);
	t36 = (-t55 * t57 + t78 * t83) * t72 * t48;
	t35 = (t50 * t55 + t61 * t83) * t48;
	t1 = [0, t36, t35, 0, 0, 0; 0, (-t59 * t72 * t38 - ((-t46 * t57 + t47 * t78) * t72 + t74 * t36) * t85) * t37, (-t53 * t38 - (t74 * t35 + t46 * t50 + t47 * t61) * t85) * t37, 0, 0, 0; 0, ((t60 * t68 + t69 * t81) * t42 - (t60 * t71 - t69 * t82) * t84) * t41, (-t42 * t71 - t68 * t84) * t54 * t41, 0, t75 * t41, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:21:55
	% EndTime: 2019-10-09 22:21:56
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (334->30), mult. (949->78), div. (65->9), fcn. (1349->13), ass. (0->50)
	t72 = sin(pkin(10));
	t74 = cos(pkin(10));
	t81 = cos(qJ(2));
	t75 = cos(pkin(6));
	t78 = sin(qJ(2));
	t85 = t75 * t78;
	t66 = t72 * t81 + t74 * t85;
	t80 = cos(qJ(3));
	t73 = sin(pkin(6));
	t77 = sin(qJ(3));
	t88 = t73 * t77;
	t59 = t66 * t80 - t74 * t88;
	t87 = t73 * t80;
	t70 = t75 * t77 + t78 * t87;
	t57 = atan2(-t59, t70);
	t54 = sin(t57);
	t55 = cos(t57);
	t48 = -t54 * t59 + t55 * t70;
	t47 = 0.1e1 / t48 ^ 2;
	t68 = -t72 * t85 + t74 * t81;
	t62 = t68 * t80 + t72 * t88;
	t93 = t47 * t62;
	t61 = t68 * t77 - t72 * t87;
	t76 = sin(qJ(5));
	t84 = t75 * t81;
	t67 = t72 * t84 + t74 * t78;
	t79 = cos(qJ(5));
	t89 = t67 * t79;
	t53 = t61 * t76 + t89;
	t51 = 0.1e1 / t53 ^ 2;
	t90 = t67 * t76;
	t52 = -t61 * t79 + t90;
	t92 = t51 * t52;
	t64 = 0.1e1 / t70 ^ 2;
	t91 = t59 * t64;
	t86 = t73 * t81;
	t83 = t52 ^ 2 * t51 + 0.1e1;
	t82 = -t54 * t70 - t55 * t59;
	t69 = t75 * t80 - t78 * t88;
	t65 = -t72 * t78 + t74 * t84;
	t63 = 0.1e1 / t70;
	t58 = t66 * t77 + t74 * t87;
	t56 = 0.1e1 / (t59 ^ 2 * t64 + 0.1e1);
	t50 = 0.1e1 / t53;
	t49 = 0.1e1 / t83;
	t46 = 0.1e1 / t48;
	t45 = 0.1e1 / (t62 ^ 2 * t47 + 0.1e1);
	t44 = (-t63 * t65 + t86 * t91) * t80 * t56;
	t43 = (t58 * t63 + t69 * t91) * t56;
	t1 = [0, t44, t43, 0, 0, 0; 0, (-t67 * t80 * t46 - ((-t54 * t65 + t55 * t86) * t80 + t82 * t44) * t93) * t45, (-t61 * t46 - (t82 * t43 + t54 * t58 + t55 * t69) * t93) * t45, 0, 0, 0; 0, ((t68 * t76 + t77 * t89) * t50 - (t68 * t79 - t77 * t90) * t92) * t49, (-t50 * t79 - t76 * t92) * t62 * t49, 0, t83 * t49, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end