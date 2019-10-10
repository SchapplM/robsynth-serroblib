% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR1
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
%   Wie in S6PRRPPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
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
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
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
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (596->31), mult. (866->78), div. (60->9), fcn. (1237->13), ass. (0->51)
	t73 = sin(pkin(10));
	t76 = cos(pkin(10));
	t79 = cos(qJ(2));
	t77 = cos(pkin(6));
	t78 = sin(qJ(2));
	t82 = t77 * t78;
	t65 = t73 * t79 + t76 * t82;
	t71 = qJ(3) + pkin(11);
	t69 = sin(t71);
	t70 = cos(t71);
	t74 = sin(pkin(6));
	t85 = t74 * t76;
	t55 = t65 * t69 + t70 * t85;
	t84 = t74 * t78;
	t62 = t69 * t84 - t77 * t70;
	t54 = atan2(-t55, t62);
	t51 = sin(t54);
	t52 = cos(t54);
	t45 = -t51 * t55 + t52 * t62;
	t44 = 0.1e1 / t45 ^ 2;
	t67 = -t73 * t82 + t76 * t79;
	t86 = t73 * t74;
	t58 = t67 * t69 - t70 * t86;
	t91 = t44 * t58;
	t59 = t67 * t70 + t69 * t86;
	t75 = cos(pkin(12));
	t81 = t77 * t79;
	t66 = t73 * t81 + t76 * t78;
	t72 = sin(pkin(12));
	t88 = t66 * t72;
	t50 = t59 * t75 + t88;
	t48 = 0.1e1 / t50 ^ 2;
	t87 = t66 * t75;
	t49 = t59 * t72 - t87;
	t90 = t48 * t49;
	t61 = 0.1e1 / t62 ^ 2;
	t89 = t55 * t61;
	t83 = t74 * t79;
	t80 = -t51 * t62 - t52 * t55;
	t64 = -t73 * t78 + t76 * t81;
	t63 = t77 * t69 + t70 * t84;
	t60 = 0.1e1 / t62;
	t57 = t65 * t70 - t69 * t85;
	t53 = 0.1e1 / (t55 ^ 2 * t61 + 0.1e1);
	t47 = 0.1e1 / t50;
	t46 = 0.1e1 / (t49 ^ 2 * t48 + 0.1e1);
	t43 = 0.1e1 / t45;
	t42 = 0.1e1 / (t58 ^ 2 * t44 + 0.1e1);
	t41 = (-t60 * t64 + t83 * t89) * t69 * t53;
	t40 = (-t57 * t60 + t63 * t89) * t53;
	t1 = [0, t41, t40, 0, 0, 0; 0, (-t66 * t69 * t43 - ((-t51 * t64 + t52 * t83) * t69 + t80 * t41) * t91) * t42, (t59 * t43 - (t80 * t40 - t51 * t57 + t52 * t63) * t91) * t42, 0, 0, 0; 0, ((-t67 * t75 - t70 * t88) * t47 - (t67 * t72 - t70 * t87) * t90) * t46, (-t47 * t72 + t75 * t90) * t58 * t46, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (689->32), mult. (949->80), div. (65->9), fcn. (1349->13), ass. (0->52)
	t80 = sin(pkin(10));
	t82 = cos(pkin(10));
	t85 = cos(qJ(2));
	t83 = cos(pkin(6));
	t84 = sin(qJ(2));
	t89 = t83 * t84;
	t70 = t80 * t85 + t82 * t89;
	t79 = qJ(3) + pkin(11);
	t75 = sin(t79);
	t77 = cos(t79);
	t81 = sin(pkin(6));
	t92 = t81 * t82;
	t60 = t70 * t75 + t77 * t92;
	t91 = t81 * t84;
	t67 = t75 * t91 - t83 * t77;
	t59 = atan2(-t60, t67);
	t56 = sin(t59);
	t57 = cos(t59);
	t50 = -t56 * t60 + t57 * t67;
	t49 = 0.1e1 / t50 ^ 2;
	t72 = -t80 * t89 + t82 * t85;
	t93 = t80 * t81;
	t63 = t72 * t75 - t77 * t93;
	t97 = t49 * t63;
	t64 = t72 * t77 + t75 * t93;
	t88 = t83 * t85;
	t71 = t80 * t88 + t82 * t84;
	t78 = pkin(12) + qJ(6);
	t74 = sin(t78);
	t76 = cos(t78);
	t55 = t64 * t76 + t71 * t74;
	t53 = 0.1e1 / t55 ^ 2;
	t54 = t64 * t74 - t71 * t76;
	t96 = t53 * t54;
	t66 = 0.1e1 / t67 ^ 2;
	t95 = t60 * t66;
	t94 = t71 * t77;
	t90 = t81 * t85;
	t87 = t54 ^ 2 * t53 + 0.1e1;
	t86 = -t56 * t67 - t57 * t60;
	t69 = -t80 * t84 + t82 * t88;
	t68 = t83 * t75 + t77 * t91;
	t65 = 0.1e1 / t67;
	t62 = t70 * t77 - t75 * t92;
	t58 = 0.1e1 / (t60 ^ 2 * t66 + 0.1e1);
	t52 = 0.1e1 / t55;
	t51 = 0.1e1 / t87;
	t48 = 0.1e1 / t50;
	t47 = 0.1e1 / (t63 ^ 2 * t49 + 0.1e1);
	t46 = (-t65 * t69 + t90 * t95) * t75 * t58;
	t45 = (-t62 * t65 + t68 * t95) * t58;
	t1 = [0, t46, t45, 0, 0, 0; 0, (-t71 * t75 * t48 - ((-t56 * t69 + t57 * t90) * t75 + t86 * t46) * t97) * t47, (t64 * t48 - (t86 * t45 - t56 * t62 + t57 * t68) * t97) * t47, 0, 0, 0; 0, ((-t72 * t76 - t74 * t94) * t52 - (t72 * t74 - t76 * t94) * t96) * t51, (-t52 * t74 + t76 * t96) * t63 * t51, 0, 0, t87 * t51;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end