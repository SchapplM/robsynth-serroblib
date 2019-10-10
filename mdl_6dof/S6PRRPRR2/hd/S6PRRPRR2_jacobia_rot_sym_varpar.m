% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR2
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
%   Wie in S6PRRPRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:28
	% EndTime: 2019-10-09 22:27:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:28
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (101->18), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->29)
	t41 = sin(pkin(11));
	t42 = sin(pkin(6));
	t53 = t41 * t42;
	t48 = cos(qJ(2));
	t52 = t42 * t48;
	t44 = cos(pkin(6));
	t46 = sin(qJ(2));
	t51 = t44 * t46;
	t50 = t44 * t48;
	t43 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (131->19), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->30)
	t50 = sin(pkin(11));
	t51 = sin(pkin(6));
	t60 = t50 * t51;
	t55 = cos(qJ(2));
	t59 = t51 * t55;
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t58 = t53 * t54;
	t57 = t53 * t55;
	t52 = cos(pkin(11));
	t43 = -t50 * t58 + t52 * t55;
	t48 = qJ(3) + pkin(12);
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
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (641->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->52)
	t71 = sin(pkin(11));
	t73 = cos(pkin(11));
	t78 = cos(qJ(2));
	t74 = cos(pkin(6));
	t76 = sin(qJ(2));
	t82 = t74 * t76;
	t64 = t71 * t78 + t73 * t82;
	t70 = qJ(3) + pkin(12);
	t68 = sin(t70);
	t69 = cos(t70);
	t72 = sin(pkin(6));
	t85 = t72 * t73;
	t54 = t64 * t68 + t69 * t85;
	t84 = t72 * t76;
	t61 = t68 * t84 - t69 * t74;
	t53 = atan2(-t54, t61);
	t48 = sin(t53);
	t49 = cos(t53);
	t44 = -t48 * t54 + t49 * t61;
	t43 = 0.1e1 / t44 ^ 2;
	t66 = -t71 * t82 + t73 * t78;
	t86 = t71 * t72;
	t57 = t66 * t68 - t69 * t86;
	t91 = t43 * t57;
	t58 = t66 * t69 + t68 * t86;
	t77 = cos(qJ(5));
	t81 = t74 * t78;
	t65 = t71 * t81 + t73 * t76;
	t75 = sin(qJ(5));
	t88 = t65 * t75;
	t51 = t58 * t77 + t88;
	t47 = 0.1e1 / t51 ^ 2;
	t87 = t65 * t77;
	t50 = t58 * t75 - t87;
	t90 = t47 * t50;
	t60 = 0.1e1 / t61 ^ 2;
	t89 = t54 * t60;
	t83 = t72 * t78;
	t80 = t47 * t50 ^ 2 + 0.1e1;
	t79 = -t48 * t61 - t49 * t54;
	t63 = -t71 * t76 + t73 * t81;
	t62 = t68 * t74 + t69 * t84;
	t59 = 0.1e1 / t61;
	t56 = t64 * t69 - t68 * t85;
	t52 = 0.1e1 / (t54 ^ 2 * t60 + 0.1e1);
	t46 = 0.1e1 / t51;
	t45 = 0.1e1 / t80;
	t42 = 0.1e1 / t44;
	t41 = 0.1e1 / (t43 * t57 ^ 2 + 0.1e1);
	t40 = (-t59 * t63 + t83 * t89) * t68 * t52;
	t39 = (-t56 * t59 + t62 * t89) * t52;
	t1 = [0, t40, t39, 0, 0, 0; 0, (-t65 * t68 * t42 - ((-t48 * t63 + t49 * t83) * t68 + t79 * t40) * t91) * t41, (t58 * t42 - (t79 * t39 - t48 * t56 + t49 * t62) * t91) * t41, 0, 0, 0; 0, ((-t66 * t77 - t69 * t88) * t46 - (t66 * t75 - t69 * t87) * t90) * t45, (-t46 * t75 + t77 * t90) * t57 * t45, 0, t80 * t45, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (748->32), mult. (1032->79), div. (70->9), fcn. (1461->13), ass. (0->54)
	t90 = sin(pkin(6));
	t91 = cos(pkin(11));
	t101 = t90 * t91;
	t89 = sin(pkin(11));
	t94 = cos(qJ(2));
	t92 = cos(pkin(6));
	t93 = sin(qJ(2));
	t98 = t92 * t93;
	t79 = t89 * t94 + t91 * t98;
	t87 = qJ(3) + pkin(12);
	t83 = sin(t87);
	t84 = cos(t87);
	t69 = t84 * t101 + t79 * t83;
	t100 = t90 * t93;
	t76 = t83 * t100 - t92 * t84;
	t68 = atan2(-t69, t76);
	t65 = sin(t68);
	t66 = cos(t68);
	t59 = -t65 * t69 + t66 * t76;
	t58 = 0.1e1 / t59 ^ 2;
	t102 = t89 * t90;
	t81 = -t89 * t98 + t91 * t94;
	t72 = -t84 * t102 + t81 * t83;
	t107 = t58 * t72;
	t97 = t92 * t94;
	t80 = t89 * t97 + t91 * t93;
	t88 = qJ(5) + qJ(6);
	t85 = sin(t88);
	t104 = t80 * t85;
	t73 = t83 * t102 + t81 * t84;
	t86 = cos(t88);
	t64 = t73 * t86 + t104;
	t62 = 0.1e1 / t64 ^ 2;
	t103 = t80 * t86;
	t63 = t73 * t85 - t103;
	t106 = t62 * t63;
	t75 = 0.1e1 / t76 ^ 2;
	t105 = t69 * t75;
	t99 = t90 * t94;
	t96 = t63 ^ 2 * t62 + 0.1e1;
	t95 = -t65 * t76 - t66 * t69;
	t78 = -t89 * t93 + t91 * t97;
	t77 = t84 * t100 + t92 * t83;
	t74 = 0.1e1 / t76;
	t71 = -t83 * t101 + t79 * t84;
	t67 = 0.1e1 / (t69 ^ 2 * t75 + 0.1e1);
	t61 = 0.1e1 / t64;
	t60 = 0.1e1 / t96;
	t57 = 0.1e1 / t59;
	t56 = 0.1e1 / (t72 ^ 2 * t58 + 0.1e1);
	t55 = (t99 * t105 - t74 * t78) * t83 * t67;
	t54 = (t77 * t105 - t71 * t74) * t67;
	t53 = t96 * t60;
	t1 = [0, t55, t54, 0, 0, 0; 0, (-t80 * t83 * t57 - ((-t65 * t78 + t66 * t99) * t83 + t95 * t55) * t107) * t56, (t73 * t57 - (t95 * t54 - t65 * t71 + t66 * t77) * t107) * t56, 0, 0, 0; 0, ((-t84 * t104 - t81 * t86) * t61 - (-t84 * t103 + t81 * t85) * t106) * t60, (t86 * t106 - t61 * t85) * t72 * t60, 0, t53, t53;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end