% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR4
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
%   Wie in S6PRRPPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPPR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (303->30), mult. (866->78), div. (60->9), fcn. (1237->13), ass. (0->48)
	t61 = sin(pkin(10));
	t64 = cos(pkin(10));
	t69 = cos(qJ(2));
	t65 = cos(pkin(6));
	t67 = sin(qJ(2));
	t72 = t65 * t67;
	t54 = t61 * t69 + t64 * t72;
	t66 = sin(qJ(3));
	t62 = sin(pkin(6));
	t68 = cos(qJ(3));
	t74 = t62 * t68;
	t46 = t54 * t66 + t64 * t74;
	t75 = t62 * t66;
	t57 = -t65 * t68 + t67 * t75;
	t45 = atan2(-t46, t57);
	t42 = sin(t45);
	t43 = cos(t45);
	t36 = -t42 * t46 + t43 * t57;
	t35 = 0.1e1 / t36 ^ 2;
	t56 = -t61 * t72 + t64 * t69;
	t49 = t56 * t66 - t61 * t74;
	t79 = t35 * t49;
	t50 = t56 * t68 + t61 * t75;
	t71 = t65 * t69;
	t55 = t61 * t71 + t64 * t67;
	t60 = sin(pkin(11));
	t63 = cos(pkin(11));
	t41 = t50 * t63 + t55 * t60;
	t39 = 0.1e1 / t41 ^ 2;
	t40 = t50 * t60 - t55 * t63;
	t78 = t39 * t40;
	t52 = 0.1e1 / t57 ^ 2;
	t77 = t46 * t52;
	t76 = t55 * t68;
	t73 = t62 * t69;
	t70 = -t42 * t57 - t43 * t46;
	t58 = t65 * t66 + t67 * t74;
	t53 = -t61 * t67 + t64 * t71;
	t51 = 0.1e1 / t57;
	t48 = t54 * t68 - t64 * t75;
	t44 = 0.1e1 / (t46 ^ 2 * t52 + 0.1e1);
	t38 = 0.1e1 / t41;
	t37 = 0.1e1 / (t40 ^ 2 * t39 + 0.1e1);
	t34 = 0.1e1 / t36;
	t33 = 0.1e1 / (t49 ^ 2 * t35 + 0.1e1);
	t32 = (-t51 * t53 + t73 * t77) * t66 * t44;
	t31 = (-t48 * t51 + t58 * t77) * t44;
	t1 = [0, t32, t31, 0, 0, 0; 0, (-t55 * t66 * t34 - ((-t42 * t53 + t43 * t73) * t66 + t70 * t32) * t79) * t33, (t50 * t34 - (t70 * t31 - t42 * t48 + t43 * t58) * t79) * t33, 0, 0, 0; 0, ((-t56 * t63 - t60 * t76) * t38 - (t56 * t60 - t63 * t76) * t78) * t37, (-t38 * t60 + t63 * t78) * t49 * t37, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (475->34), mult. (1334->86), div. (60->9), fcn. (1876->13), ass. (0->50)
	t77 = sin(pkin(10));
	t80 = cos(pkin(10));
	t85 = cos(qJ(2));
	t81 = cos(pkin(6));
	t83 = sin(qJ(2));
	t89 = t81 * t83;
	t72 = t77 * t85 + t80 * t89;
	t76 = sin(pkin(11));
	t79 = cos(pkin(11));
	t84 = cos(qJ(3));
	t88 = t81 * t85;
	t86 = -t77 * t83 + t80 * t88;
	t78 = sin(pkin(6));
	t82 = sin(qJ(3));
	t91 = t78 * t82;
	t58 = (t72 * t84 - t80 * t91) * t76 + t86 * t79;
	t90 = t78 * t84;
	t66 = (t81 * t82 + t83 * t90) * t76 + t78 * t85 * t79;
	t55 = atan2(-t58, t66);
	t51 = sin(t55);
	t52 = cos(t55);
	t50 = -t51 * t58 + t52 * t66;
	t49 = 0.1e1 / t50 ^ 2;
	t74 = -t77 * t89 + t80 * t85;
	t69 = t74 * t84 + t77 * t91;
	t73 = t77 * t88 + t80 * t83;
	t93 = t73 * t79;
	t60 = t69 * t76 - t93;
	t96 = t49 * t60;
	t64 = 0.1e1 / t66 ^ 2;
	t95 = t58 * t64;
	t61 = t69 * t79 + t73 * t76;
	t57 = 0.1e1 / t61 ^ 2;
	t68 = -t74 * t82 + t77 * t90;
	t94 = t68 ^ 2 * t57;
	t92 = t76 * t84;
	t87 = -t51 * t66 - t52 * t58;
	t75 = t81 * t84 - t83 * t91;
	t70 = (-t79 * t83 + t85 * t92) * t78;
	t67 = -t72 * t82 - t80 * t90;
	t63 = 0.1e1 / t66;
	t62 = -t72 * t79 + t86 * t92;
	t56 = 0.1e1 / t61;
	t54 = 0.1e1 / (t58 ^ 2 * t64 + 0.1e1);
	t53 = 0.1e1 / (0.1e1 + t94);
	t48 = 0.1e1 / t50;
	t47 = 0.1e1 / (t60 ^ 2 * t49 + 0.1e1);
	t46 = (-t63 * t67 + t75 * t95) * t76 * t54;
	t45 = (-t62 * t63 + t70 * t95) * t54;
	t1 = [0, t45, t46, 0, 0, 0; 0, ((-t73 * t92 - t74 * t79) * t48 - (t87 * t45 - t51 * t62 + t52 * t70) * t96) * t47, (t68 * t76 * t48 - ((-t51 * t67 + t52 * t75) * t76 + t87 * t46) * t96) * t47, 0, 0, 0; 0, (t73 * t82 * t56 - (t74 * t76 - t84 * t93) * t68 * t57) * t53, (-t56 * t69 - t79 * t94) * t53, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:44
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (453->35), mult. (1262->91), div. (65->9), fcn. (1781->15), ass. (0->55)
	t88 = sin(pkin(6));
	t96 = cos(qJ(3));
	t103 = t88 * t96;
	t91 = cos(pkin(6));
	t94 = sin(qJ(2));
	t101 = t91 * t94;
	t87 = sin(pkin(10));
	t90 = cos(pkin(10));
	t97 = cos(qJ(2));
	t81 = t90 * t101 + t87 * t97;
	t93 = sin(qJ(3));
	t74 = t90 * t103 + t81 * t93;
	t104 = t88 * t93;
	t84 = -t94 * t104 + t91 * t96;
	t71 = atan2(t74, t84);
	t68 = sin(t71);
	t69 = cos(t71);
	t61 = t68 * t74 + t69 * t84;
	t60 = 0.1e1 / t61 ^ 2;
	t83 = -t87 * t101 + t90 * t97;
	t76 = t87 * t103 - t83 * t93;
	t108 = t60 * t76;
	t77 = t87 * t104 + t83 * t96;
	t100 = t91 * t97;
	t82 = t87 * t100 + t90 * t94;
	t86 = sin(pkin(11));
	t89 = cos(pkin(11));
	t66 = t77 * t86 - t82 * t89;
	t67 = t77 * t89 + t82 * t86;
	t92 = sin(qJ(6));
	t95 = cos(qJ(6));
	t65 = t66 * t92 + t67 * t95;
	t63 = 0.1e1 / t65 ^ 2;
	t64 = -t66 * t95 + t67 * t92;
	t107 = t63 * t64;
	t79 = 0.1e1 / t84 ^ 2;
	t106 = t74 * t79;
	t105 = t82 * t96;
	t102 = t88 * t97;
	t99 = t64 ^ 2 * t63 + 0.1e1;
	t98 = -t68 * t84 + t69 * t74;
	t85 = -t94 * t103 - t91 * t93;
	t80 = t90 * t100 - t87 * t94;
	t78 = 0.1e1 / t84;
	t75 = -t90 * t104 + t81 * t96;
	t73 = -t89 * t105 + t83 * t86;
	t72 = -t86 * t105 - t83 * t89;
	t70 = 0.1e1 / (t74 ^ 2 * t79 + 0.1e1);
	t62 = 0.1e1 / t65;
	t59 = 0.1e1 / t61;
	t58 = 0.1e1 / (t76 ^ 2 * t60 + 0.1e1);
	t57 = (t102 * t106 + t78 * t80) * t93 * t70;
	t56 = (-t85 * t106 + t75 * t78) * t70;
	t55 = 0.1e1 / t99;
	t1 = [0, t57, t56, 0, 0, 0; 0, (t82 * t93 * t59 - ((-t69 * t102 + t68 * t80) * t93 + t98 * t57) * t108) * t58, (-t77 * t59 - (t98 * t56 + t68 * t75 + t69 * t85) * t108) * t58, 0, 0, 0; 0, ((-t72 * t95 + t73 * t92) * t62 - (t72 * t92 + t73 * t95) * t107) * t55, ((-t86 * t95 + t89 * t92) * t62 - (t86 * t92 + t89 * t95) * t107) * t55 * t76, 0, 0, t99 * t55;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end