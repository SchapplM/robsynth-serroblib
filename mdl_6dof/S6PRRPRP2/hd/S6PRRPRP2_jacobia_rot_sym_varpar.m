% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRP2
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
%   Wie in S6PRRPRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:18
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRP2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRP2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:18:13
	% EndTime: 2019-10-09 22:18:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:14
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
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:14
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
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:14
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (641->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->52)
	t71 = sin(pkin(10));
	t73 = cos(pkin(10));
	t78 = cos(qJ(2));
	t74 = cos(pkin(6));
	t76 = sin(qJ(2));
	t82 = t74 * t76;
	t64 = t71 * t78 + t73 * t82;
	t70 = qJ(3) + pkin(11);
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
	% StartTime: 2019-10-09 22:18:14
	% EndTime: 2019-10-09 22:18:14
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (1149->41), mult. (1993->102), div. (87->9), fcn. (2807->13), ass. (0->60)
	t101 = sin(qJ(5));
	t102 = sin(qJ(2));
	t100 = cos(pkin(6));
	t104 = cos(qJ(2));
	t108 = t100 * t104;
	t97 = sin(pkin(10));
	t99 = cos(pkin(10));
	t105 = -t97 * t102 + t99 * t108;
	t119 = t105 * t101;
	t103 = cos(qJ(5));
	t98 = sin(pkin(6));
	t113 = t98 * t99;
	t109 = t100 * t102;
	t91 = t97 * t104 + t99 * t109;
	t96 = qJ(3) + pkin(11);
	t94 = sin(t96);
	t95 = cos(t96);
	t80 = -t94 * t113 + t91 * t95;
	t72 = t80 * t101 + t105 * t103;
	t112 = t102 * t98;
	t89 = t100 * t94 + t95 * t112;
	t85 = t98 * t104 * t103 + t89 * t101;
	t69 = atan2(-t72, t85);
	t66 = sin(t69);
	t67 = cos(t69);
	t64 = -t66 * t72 + t67 * t85;
	t63 = 0.1e1 / t64 ^ 2;
	t92 = t99 * t102 + t97 * t108;
	t110 = t92 * t103;
	t114 = t97 * t98;
	t93 = t99 * t104 - t97 * t109;
	t82 = t94 * t114 + t93 * t95;
	t75 = t82 * t101 - t110;
	t118 = t63 * t75;
	t111 = t92 * t101;
	t76 = t82 * t103 + t111;
	t71 = 0.1e1 / t76 ^ 2;
	t81 = t95 * t114 - t93 * t94;
	t117 = t71 * t81;
	t84 = 0.1e1 / t85 ^ 2;
	t116 = t72 * t84;
	t115 = t81 ^ 2 * t71;
	t107 = t101 * t104;
	t106 = -t66 * t85 - t67 * t72;
	t88 = t100 * t95 - t94 * t112;
	t87 = (-t102 * t103 + t95 * t107) * t98;
	t86 = t89 * t103 - t98 * t107;
	t83 = 0.1e1 / t85;
	t79 = -t95 * t113 - t91 * t94;
	t77 = -t91 * t103 + t95 * t119;
	t74 = t80 * t103 - t119;
	t70 = 0.1e1 / t76;
	t68 = 0.1e1 / (t72 ^ 2 * t84 + 0.1e1);
	t65 = 0.1e1 / (0.1e1 + t115);
	t62 = 0.1e1 / t64;
	t61 = 0.1e1 / (t75 ^ 2 * t63 + 0.1e1);
	t60 = (t88 * t116 - t79 * t83) * t68 * t101;
	t59 = (t87 * t116 - t77 * t83) * t68;
	t58 = (t86 * t116 - t74 * t83) * t68;
	t1 = [0, t59, t60, 0, t58, 0; 0, ((-t93 * t103 - t95 * t111) * t62 - (t106 * t59 - t66 * t77 + t67 * t87) * t118) * t61, (t81 * t101 * t62 - (t106 * t60 + (-t66 * t79 + t67 * t88) * t101) * t118) * t61, 0, (t76 * t62 - (t106 * t58 - t66 * t74 + t67 * t86) * t118) * t61, 0; 0, (t92 * t94 * t70 - (t93 * t101 - t95 * t110) * t117) * t65, (-t103 * t115 - t70 * t82) * t65, 0, t75 * t65 * t117, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end