% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPP1
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
%   Wie in S6PRRRPP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
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
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (334->30), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->49)
	t62 = sin(pkin(10));
	t64 = cos(pkin(10));
	t71 = cos(qJ(2));
	t65 = cos(pkin(6));
	t68 = sin(qJ(2));
	t75 = t65 * t68;
	t56 = t62 * t71 + t64 * t75;
	t67 = sin(qJ(3));
	t63 = sin(pkin(6));
	t70 = cos(qJ(3));
	t77 = t63 * t70;
	t48 = t56 * t67 + t64 * t77;
	t78 = t63 * t67;
	t59 = -t65 * t70 + t68 * t78;
	t47 = atan2(-t48, t59);
	t44 = sin(t47);
	t45 = cos(t47);
	t38 = -t44 * t48 + t45 * t59;
	t37 = 0.1e1 / t38 ^ 2;
	t58 = -t62 * t75 + t64 * t71;
	t51 = t58 * t67 - t62 * t77;
	t82 = t37 * t51;
	t52 = t58 * t70 + t62 * t78;
	t74 = t65 * t71;
	t57 = t62 * t74 + t64 * t68;
	t66 = sin(qJ(4));
	t69 = cos(qJ(4));
	t43 = t52 * t69 + t57 * t66;
	t41 = 0.1e1 / t43 ^ 2;
	t42 = t52 * t66 - t57 * t69;
	t81 = t41 * t42;
	t54 = 0.1e1 / t59 ^ 2;
	t80 = t48 * t54;
	t79 = t57 * t70;
	t76 = t63 * t71;
	t73 = t42 ^ 2 * t41 + 0.1e1;
	t72 = -t44 * t59 - t45 * t48;
	t60 = t65 * t67 + t68 * t77;
	t55 = -t62 * t68 + t64 * t74;
	t53 = 0.1e1 / t59;
	t50 = t56 * t70 - t64 * t78;
	t46 = 0.1e1 / (t48 ^ 2 * t54 + 0.1e1);
	t40 = 0.1e1 / t43;
	t39 = 0.1e1 / t73;
	t36 = 0.1e1 / t38;
	t35 = 0.1e1 / (t51 ^ 2 * t37 + 0.1e1);
	t34 = (-t53 * t55 + t76 * t80) * t67 * t46;
	t33 = (-t50 * t53 + t60 * t80) * t46;
	t1 = [0, t34, t33, 0, 0, 0; 0, (-t57 * t67 * t36 - ((-t44 * t55 + t45 * t76) * t67 + t72 * t34) * t82) * t35, (t52 * t36 - (t72 * t33 - t44 * t50 + t45 * t60) * t82) * t35, 0, 0, 0; 0, ((-t58 * t69 - t66 * t79) * t40 - (t58 * t66 - t69 * t79) * t81) * t39, (-t40 * t66 + t69 * t81) * t51 * t39, t73 * t39, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (382->31), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->50)
	t72 = sin(pkin(10));
	t74 = cos(pkin(10));
	t79 = cos(qJ(2));
	t75 = cos(pkin(6));
	t77 = sin(qJ(2));
	t83 = t75 * t77;
	t63 = t72 * t79 + t74 * t83;
	t76 = sin(qJ(3));
	t73 = sin(pkin(6));
	t78 = cos(qJ(3));
	t85 = t73 * t78;
	t55 = t63 * t76 + t74 * t85;
	t86 = t73 * t76;
	t66 = -t75 * t78 + t77 * t86;
	t54 = atan2(-t55, t66);
	t51 = sin(t54);
	t52 = cos(t54);
	t45 = -t51 * t55 + t52 * t66;
	t44 = 0.1e1 / t45 ^ 2;
	t65 = -t72 * t83 + t74 * t79;
	t58 = t65 * t76 - t72 * t85;
	t90 = t44 * t58;
	t59 = t65 * t78 + t72 * t86;
	t82 = t75 * t79;
	t64 = t72 * t82 + t74 * t77;
	t71 = qJ(4) + pkin(11);
	t69 = sin(t71);
	t70 = cos(t71);
	t50 = t59 * t70 + t64 * t69;
	t48 = 0.1e1 / t50 ^ 2;
	t49 = t59 * t69 - t64 * t70;
	t89 = t48 * t49;
	t61 = 0.1e1 / t66 ^ 2;
	t88 = t55 * t61;
	t87 = t64 * t78;
	t84 = t73 * t79;
	t81 = t49 ^ 2 * t48 + 0.1e1;
	t80 = -t51 * t66 - t52 * t55;
	t67 = t75 * t76 + t77 * t85;
	t62 = -t72 * t77 + t74 * t82;
	t60 = 0.1e1 / t66;
	t57 = t63 * t78 - t74 * t86;
	t53 = 0.1e1 / (t55 ^ 2 * t61 + 0.1e1);
	t47 = 0.1e1 / t50;
	t46 = 0.1e1 / t81;
	t43 = 0.1e1 / t45;
	t42 = 0.1e1 / (t58 ^ 2 * t44 + 0.1e1);
	t41 = (-t60 * t62 + t84 * t88) * t76 * t53;
	t40 = (-t57 * t60 + t67 * t88) * t53;
	t1 = [0, t41, t40, 0, 0, 0; 0, (-t64 * t76 * t43 - ((-t51 * t62 + t52 * t84) * t76 + t80 * t41) * t90) * t42, (t59 * t43 - (t80 * t40 - t51 * t57 + t52 * t67) * t90) * t42, 0, 0, 0; 0, ((-t65 * t70 - t69 * t87) * t47 - (t65 * t69 - t70 * t87) * t89) * t46, (-t47 * t69 + t70 * t89) * t58 * t46, t81 * t46, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (1132->41), mult. (1993->101), div. (87->9), fcn. (2807->13), ass. (0->58)
	t101 = cos(pkin(10));
	t104 = sin(qJ(2));
	t102 = cos(pkin(6));
	t106 = cos(qJ(2));
	t109 = t102 * t106;
	t99 = sin(pkin(10));
	t107 = t101 * t109 - t99 * t104;
	t105 = cos(qJ(3));
	t100 = sin(pkin(6));
	t103 = sin(qJ(3));
	t113 = t100 * t103;
	t110 = t102 * t104;
	t91 = t101 * t110 + t99 * t106;
	t86 = -t101 * t113 + t91 * t105;
	t98 = qJ(4) + pkin(11);
	t96 = sin(t98);
	t97 = cos(t98);
	t74 = t107 * t97 + t86 * t96;
	t111 = t100 * t106;
	t112 = t100 * t105;
	t95 = t102 * t103 + t104 * t112;
	t82 = t97 * t111 + t95 * t96;
	t70 = atan2(-t74, t82);
	t67 = sin(t70);
	t68 = cos(t70);
	t66 = -t67 * t74 + t68 * t82;
	t65 = 0.1e1 / t66 ^ 2;
	t92 = t101 * t104 + t99 * t109;
	t115 = t92 * t97;
	t93 = t101 * t106 - t99 * t110;
	t88 = t93 * t105 + t99 * t113;
	t77 = t88 * t96 - t115;
	t119 = t65 * t77;
	t78 = t88 * t97 + t92 * t96;
	t73 = 0.1e1 / t78 ^ 2;
	t87 = -t93 * t103 + t99 * t112;
	t118 = t73 * t87;
	t81 = 0.1e1 / t82 ^ 2;
	t117 = t74 * t81;
	t116 = t87 ^ 2 * t73;
	t114 = t105 * t96;
	t108 = -t67 * t82 - t68 * t74;
	t94 = t102 * t105 - t104 * t113;
	t89 = (-t104 * t97 + t106 * t114) * t100;
	t85 = -t101 * t112 - t91 * t103;
	t83 = -t96 * t111 + t95 * t97;
	t80 = 0.1e1 / t82;
	t79 = t107 * t114 - t91 * t97;
	t76 = -t107 * t96 + t86 * t97;
	t72 = 0.1e1 / t78;
	t71 = 0.1e1 / (0.1e1 + t116);
	t69 = 0.1e1 / (t74 ^ 2 * t81 + 0.1e1);
	t64 = 0.1e1 / t66;
	t63 = 0.1e1 / (t77 ^ 2 * t65 + 0.1e1);
	t62 = (t94 * t117 - t80 * t85) * t96 * t69;
	t61 = (t89 * t117 - t79 * t80) * t69;
	t60 = (t83 * t117 - t76 * t80) * t69;
	t1 = [0, t61, t62, t60, 0, 0; 0, ((-t92 * t114 - t93 * t97) * t64 - (t108 * t61 - t67 * t79 + t68 * t89) * t119) * t63, (t87 * t96 * t64 - ((-t67 * t85 + t68 * t94) * t96 + t108 * t62) * t119) * t63, (t78 * t64 - (t108 * t60 - t67 * t76 + t68 * t83) * t119) * t63, 0, 0; 0, (t92 * t103 * t72 - (-t105 * t115 + t93 * t96) * t118) * t71, (-t97 * t116 - t72 * t88) * t71, t77 * t71 * t118, 0, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end