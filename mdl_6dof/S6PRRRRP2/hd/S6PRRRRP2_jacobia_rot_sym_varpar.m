% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRP2
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
%   Wie in S6PRRRRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
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
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
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
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:56
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
	% StartTime: 2019-10-09 23:03:56
	% EndTime: 2019-10-09 23:03:57
	% DurationCPUTime: 0.36s
	% Computational Cost: add. (1533->43), mult. (2650->105), div. (117->9), fcn. (3731->13), ass. (0->60)
	t130 = sin(pkin(11));
	t132 = cos(pkin(11));
	t137 = cos(qJ(2));
	t133 = cos(pkin(6));
	t135 = sin(qJ(2));
	t141 = t133 * t135;
	t126 = -t130 * t141 + t132 * t137;
	t129 = qJ(3) + qJ(4);
	t127 = sin(t129);
	t128 = cos(t129);
	t131 = sin(pkin(6));
	t144 = t130 * t131;
	t115 = t126 * t128 + t127 * t144;
	t134 = sin(qJ(5));
	t140 = t133 * t137;
	t125 = t130 * t140 + t132 * t135;
	t136 = cos(qJ(5));
	t146 = t125 * t136;
	t108 = t115 * t134 - t146;
	t124 = t130 * t137 + t132 * t141;
	t143 = t131 * t132;
	t113 = t124 * t128 - t127 * t143;
	t138 = -t130 * t135 + t132 * t140;
	t105 = t113 * t134 + t138 * t136;
	t142 = t131 * t135;
	t122 = t133 * t127 + t128 * t142;
	t118 = t131 * t137 * t136 + t122 * t134;
	t102 = atan2(-t105, t118);
	t100 = cos(t102);
	t99 = sin(t102);
	t97 = t100 * t118 - t99 * t105;
	t96 = 0.1e1 / t97 ^ 2;
	t150 = t108 * t96;
	t109 = t115 * t136 + t125 * t134;
	t104 = 0.1e1 / t109 ^ 2;
	t114 = -t126 * t127 + t128 * t144;
	t149 = t104 * t114;
	t117 = 0.1e1 / t118 ^ 2;
	t148 = t105 * t117;
	t147 = t114 ^ 2 * t104;
	t145 = t128 * t134;
	t139 = t134 * t137;
	t121 = -t127 * t142 + t133 * t128;
	t120 = (t128 * t139 - t135 * t136) * t131;
	t119 = t122 * t136 - t131 * t139;
	t116 = 0.1e1 / t118;
	t112 = -t124 * t127 - t128 * t143;
	t110 = -t124 * t136 + t138 * t145;
	t107 = t113 * t136 - t138 * t134;
	t103 = 0.1e1 / t109;
	t101 = 0.1e1 / (t105 ^ 2 * t117 + 0.1e1);
	t98 = 0.1e1 / (0.1e1 + t147);
	t95 = 0.1e1 / t97;
	t94 = 0.1e1 / (t108 ^ 2 * t96 + 0.1e1);
	t93 = (-t112 * t116 + t121 * t148) * t134 * t101;
	t92 = (-t103 * t115 - t136 * t147) * t98;
	t91 = (-t110 * t116 + t120 * t148) * t101;
	t90 = (-t107 * t116 + t119 * t148) * t101;
	t89 = (t114 * t134 * t95 - ((-t112 * t134 - t118 * t93) * t99 + (-t105 * t93 + t121 * t134) * t100) * t150) * t94;
	t1 = [0, t91, t93, t93, t90, 0; 0, ((-t125 * t145 - t126 * t136) * t95 - ((-t118 * t91 - t110) * t99 + (-t105 * t91 + t120) * t100) * t150) * t94, t89, t89, (t109 * t95 - ((-t118 * t90 - t107) * t99 + (-t105 * t90 + t119) * t100) * t150) * t94, 0; 0, (t125 * t127 * t103 - (t126 * t134 - t128 * t146) * t149) * t98, t92, t92, t108 * t98 * t149, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end