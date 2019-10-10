% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRP1
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
%   Wie in S6PRRRRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (162->19), mult. (316->43), div. (52->11), fcn. (458->11), ass. (0->31)
	t59 = sin(pkin(11));
	t60 = sin(pkin(6));
	t69 = t59 * t60;
	t64 = cos(qJ(2));
	t68 = t60 * t64;
	t62 = cos(pkin(6));
	t63 = sin(qJ(2));
	t67 = t62 * t63;
	t66 = t62 * t64;
	t61 = cos(pkin(11));
	t52 = -t59 * t67 + t61 * t64;
	t58 = qJ(3) + qJ(4);
	t54 = sin(t58);
	t55 = cos(t58);
	t43 = t52 * t55 + t54 * t69;
	t41 = 0.1e1 / t43 ^ 2;
	t42 = t52 * t54 - t55 * t69;
	t65 = t42 ^ 2 * t41 + 0.1e1;
	t57 = 0.1e1 / t64 ^ 2;
	t51 = t59 * t66 + t61 * t63;
	t50 = t59 * t64 + t61 * t67;
	t48 = t59 * t63 - t61 * t66;
	t46 = atan2(-t48, -t68);
	t45 = cos(t46);
	t44 = sin(t46);
	t40 = 0.1e1 / t65;
	t39 = -t44 * t48 - t45 * t68;
	t38 = 0.1e1 / t39 ^ 2;
	t36 = (t50 / t64 + t63 * t48 * t57) / t60 / (0.1e1 + t48 ^ 2 / t60 ^ 2 * t57);
	t35 = t65 * t40;
	t1 = [0, t36, 0, 0, 0, 0; 0, (t52 / t39 - (t45 * t60 * t63 - t44 * t50 + (t44 * t68 - t45 * t48) * t36) * t51 * t38) / (t51 ^ 2 * t38 + 0.1e1), 0, 0, 0, 0; 0, (-t54 / t43 + t55 * t42 * t41) * t51 * t40, t35, t35, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (948->31), mult. (1387->79), div. (95->9), fcn. (1976->13), ass. (0->54)
	t92 = sin(pkin(6));
	t93 = cos(pkin(11));
	t105 = t92 * t93;
	t94 = cos(pkin(6));
	t96 = sin(qJ(2));
	t102 = t94 * t96;
	t91 = sin(pkin(11));
	t98 = cos(qJ(2));
	t84 = t93 * t102 + t91 * t98;
	t90 = qJ(3) + qJ(4);
	t88 = sin(t90);
	t89 = cos(t90);
	t74 = t89 * t105 + t84 * t88;
	t104 = t92 * t96;
	t81 = t88 * t104 - t94 * t89;
	t73 = atan2(-t74, t81);
	t70 = sin(t73);
	t71 = cos(t73);
	t64 = -t70 * t74 + t71 * t81;
	t63 = 0.1e1 / t64 ^ 2;
	t106 = t91 * t92;
	t86 = -t91 * t102 + t93 * t98;
	t77 = -t89 * t106 + t86 * t88;
	t111 = t63 * t77;
	t101 = t94 * t98;
	t85 = t91 * t101 + t93 * t96;
	t95 = sin(qJ(5));
	t108 = t85 * t95;
	t78 = t88 * t106 + t86 * t89;
	t97 = cos(qJ(5));
	t69 = t78 * t97 + t108;
	t67 = 0.1e1 / t69 ^ 2;
	t107 = t85 * t97;
	t68 = t78 * t95 - t107;
	t110 = t67 * t68;
	t80 = 0.1e1 / t81 ^ 2;
	t109 = t74 * t80;
	t103 = t92 * t98;
	t100 = t68 ^ 2 * t67 + 0.1e1;
	t99 = -t70 * t81 - t71 * t74;
	t83 = t93 * t101 - t91 * t96;
	t82 = t89 * t104 + t94 * t88;
	t79 = 0.1e1 / t81;
	t76 = -t88 * t105 + t84 * t89;
	t72 = 0.1e1 / (t74 ^ 2 * t80 + 0.1e1);
	t66 = 0.1e1 / t69;
	t65 = 0.1e1 / t100;
	t62 = 0.1e1 / t64;
	t61 = 0.1e1 / (t77 ^ 2 * t63 + 0.1e1);
	t60 = (t103 * t109 - t79 * t83) * t88 * t72;
	t59 = (t82 * t109 - t76 * t79) * t72;
	t58 = (t97 * t110 - t66 * t95) * t77 * t65;
	t57 = (t78 * t62 - (t99 * t59 - t70 * t76 + t71 * t82) * t111) * t61;
	t1 = [0, t60, t59, t59, 0, 0; 0, (-t85 * t88 * t62 - ((t71 * t103 - t70 * t83) * t88 + t99 * t60) * t111) * t61, t57, t57, 0, 0; 0, ((-t89 * t108 - t86 * t97) * t66 - (-t89 * t107 + t86 * t95) * t110) * t65, t58, t58, t100 * t65, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:02:05
	% EndTime: 2019-10-09 23:02:05
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (948->31), mult. (1387->79), div. (95->9), fcn. (1976->13), ass. (0->54)
	t101 = sin(pkin(6));
	t102 = cos(pkin(11));
	t114 = t101 * t102;
	t100 = sin(pkin(11));
	t107 = cos(qJ(2));
	t103 = cos(pkin(6));
	t105 = sin(qJ(2));
	t111 = t103 * t105;
	t93 = t100 * t107 + t102 * t111;
	t99 = qJ(3) + qJ(4);
	t97 = sin(t99);
	t98 = cos(t99);
	t83 = t98 * t114 + t93 * t97;
	t113 = t101 * t105;
	t90 = -t103 * t98 + t97 * t113;
	t82 = atan2(-t83, t90);
	t79 = sin(t82);
	t80 = cos(t82);
	t73 = -t79 * t83 + t80 * t90;
	t72 = 0.1e1 / t73 ^ 2;
	t115 = t100 * t101;
	t95 = -t100 * t111 + t102 * t107;
	t86 = -t98 * t115 + t95 * t97;
	t120 = t72 * t86;
	t106 = cos(qJ(5));
	t104 = sin(qJ(5));
	t110 = t103 * t107;
	t94 = t100 * t110 + t102 * t105;
	t117 = t94 * t104;
	t87 = t97 * t115 + t95 * t98;
	t78 = t87 * t106 + t117;
	t76 = 0.1e1 / t78 ^ 2;
	t116 = t94 * t106;
	t77 = t87 * t104 - t116;
	t119 = t76 * t77;
	t89 = 0.1e1 / t90 ^ 2;
	t118 = t83 * t89;
	t112 = t101 * t107;
	t109 = t77 ^ 2 * t76 + 0.1e1;
	t108 = -t79 * t90 - t80 * t83;
	t92 = -t100 * t105 + t102 * t110;
	t91 = t103 * t97 + t98 * t113;
	t88 = 0.1e1 / t90;
	t85 = -t97 * t114 + t93 * t98;
	t81 = 0.1e1 / (t83 ^ 2 * t89 + 0.1e1);
	t75 = 0.1e1 / t78;
	t74 = 0.1e1 / t109;
	t71 = 0.1e1 / t73;
	t70 = 0.1e1 / (t86 ^ 2 * t72 + 0.1e1);
	t69 = (t112 * t118 - t88 * t92) * t97 * t81;
	t68 = (t91 * t118 - t85 * t88) * t81;
	t67 = (-t104 * t75 + t106 * t119) * t86 * t74;
	t66 = (t87 * t71 - (t108 * t68 - t79 * t85 + t80 * t91) * t120) * t70;
	t1 = [0, t69, t68, t68, 0, 0; 0, (-t94 * t97 * t71 - ((t80 * t112 - t79 * t92) * t97 + t108 * t69) * t120) * t70, t66, t66, 0, 0; 0, ((-t95 * t106 - t98 * t117) * t75 - (t95 * t104 - t98 * t116) * t119) * t74, t67, t67, t109 * t74, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end