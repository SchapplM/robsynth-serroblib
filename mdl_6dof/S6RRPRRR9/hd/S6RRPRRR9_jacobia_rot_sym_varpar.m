% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR9
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
%   Wie in S6RRPRRR9_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:03
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR9_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR9_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
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
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (156->24), mult. (403->63), div. (67->11), fcn. (610->11), ass. (0->39)
	t46 = sin(qJ(2));
	t48 = cos(qJ(2));
	t49 = cos(qJ(1));
	t47 = sin(qJ(1));
	t52 = cos(pkin(6));
	t51 = t47 * t52;
	t38 = -t46 * t51 + t49 * t48;
	t43 = sin(pkin(12));
	t45 = cos(pkin(12));
	t44 = sin(pkin(6));
	t55 = t44 * t47;
	t29 = t38 * t45 + t43 * t55;
	t27 = 0.1e1 / t29 ^ 2;
	t28 = t38 * t43 - t45 * t55;
	t59 = t27 * t28;
	t50 = t49 * t52;
	t34 = t47 * t46 - t48 * t50;
	t54 = t44 * t48;
	t32 = atan2(-t34, -t54);
	t31 = cos(t32);
	t58 = t31 * t34;
	t30 = sin(t32);
	t24 = -t30 * t34 - t31 * t54;
	t23 = 0.1e1 / t24 ^ 2;
	t37 = t49 * t46 + t48 * t51;
	t57 = t37 ^ 2 * t23;
	t40 = 0.1e1 / t44;
	t41 = 0.1e1 / t48;
	t56 = t40 * t41;
	t53 = t44 * t49;
	t42 = 0.1e1 / t48 ^ 2;
	t36 = t46 * t50 + t47 * t48;
	t33 = 0.1e1 / (0.1e1 + t34 ^ 2 / t44 ^ 2 * t42);
	t26 = 0.1e1 / t29;
	t25 = 0.1e1 / (t28 ^ 2 * t27 + 0.1e1);
	t22 = 0.1e1 / t24;
	t21 = 0.1e1 / (0.1e1 + t57);
	t20 = (t34 * t42 * t46 + t36 * t41) * t40 * t33;
	t1 = [t37 * t33 * t56, t20, 0, 0, 0, 0; (-t34 * t22 - (-t30 + (-t56 * t58 + t30) * t33) * t57) * t21, (t38 * t22 - (t31 * t44 * t46 - t30 * t36 + (t30 * t54 - t58) * t20) * t37 * t23) * t21, 0, 0, 0, 0; ((-t36 * t43 - t45 * t53) * t26 - (-t36 * t45 + t43 * t53) * t59) * t25, (-t43 * t26 + t45 * t59) * t37 * t25, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (221->25), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->41)
	t58 = sin(qJ(2));
	t60 = cos(qJ(2));
	t61 = cos(qJ(1));
	t59 = sin(qJ(1));
	t65 = cos(pkin(6));
	t63 = t59 * t65;
	t49 = -t58 * t63 + t61 * t60;
	t54 = pkin(12) + qJ(4);
	t51 = sin(t54);
	t52 = cos(t54);
	t57 = sin(pkin(6));
	t68 = t57 * t59;
	t40 = t49 * t52 + t51 * t68;
	t38 = 0.1e1 / t40 ^ 2;
	t39 = t49 * t51 - t52 * t68;
	t72 = t38 * t39;
	t62 = t61 * t65;
	t45 = t59 * t58 - t60 * t62;
	t67 = t57 * t60;
	t43 = atan2(-t45, -t67);
	t42 = cos(t43);
	t71 = t42 * t45;
	t41 = sin(t43);
	t35 = -t41 * t45 - t42 * t67;
	t34 = 0.1e1 / t35 ^ 2;
	t48 = t61 * t58 + t60 * t63;
	t70 = t48 ^ 2 * t34;
	t53 = 0.1e1 / t57;
	t55 = 0.1e1 / t60;
	t69 = t53 * t55;
	t66 = t57 * t61;
	t64 = t39 ^ 2 * t38 + 0.1e1;
	t56 = 0.1e1 / t60 ^ 2;
	t47 = t58 * t62 + t59 * t60;
	t44 = 0.1e1 / (0.1e1 + t45 ^ 2 / t57 ^ 2 * t56);
	t37 = 0.1e1 / t40;
	t36 = 0.1e1 / t64;
	t33 = 0.1e1 / t35;
	t32 = 0.1e1 / (0.1e1 + t70);
	t31 = (t45 * t56 * t58 + t47 * t55) * t53 * t44;
	t1 = [t48 * t44 * t69, t31, 0, 0, 0, 0; (-t45 * t33 - (-t41 + (-t69 * t71 + t41) * t44) * t70) * t32, (t49 * t33 - (t42 * t57 * t58 - t41 * t47 + (t41 * t67 - t71) * t31) * t48 * t34) * t32, 0, 0, 0, 0; ((-t47 * t51 - t52 * t66) * t37 - (-t47 * t52 + t51 * t66) * t72) * t36, (-t51 * t37 + t52 * t72) * t48 * t36, 0, t64 * t36, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (314->25), mult. (499->64), div. (77->11), fcn. (736->11), ass. (0->42)
	t71 = sin(qJ(2));
	t73 = cos(qJ(2));
	t74 = cos(qJ(1));
	t72 = sin(qJ(1));
	t78 = cos(pkin(6));
	t76 = t72 * t78;
	t62 = -t71 * t76 + t74 * t73;
	t66 = pkin(12) + qJ(4) + qJ(5);
	t64 = sin(t66);
	t65 = cos(t66);
	t70 = sin(pkin(6));
	t81 = t70 * t72;
	t53 = t62 * t65 + t64 * t81;
	t51 = 0.1e1 / t53 ^ 2;
	t52 = t62 * t64 - t65 * t81;
	t85 = t51 * t52;
	t75 = t74 * t78;
	t58 = t72 * t71 - t73 * t75;
	t80 = t70 * t73;
	t56 = atan2(-t58, -t80);
	t55 = cos(t56);
	t84 = t55 * t58;
	t54 = sin(t56);
	t48 = -t54 * t58 - t55 * t80;
	t47 = 0.1e1 / t48 ^ 2;
	t61 = t74 * t71 + t73 * t76;
	t83 = t61 ^ 2 * t47;
	t67 = 0.1e1 / t70;
	t68 = 0.1e1 / t73;
	t82 = t67 * t68;
	t79 = t70 * t74;
	t77 = t52 ^ 2 * t51 + 0.1e1;
	t69 = 0.1e1 / t73 ^ 2;
	t60 = t71 * t75 + t72 * t73;
	t57 = 0.1e1 / (0.1e1 + t58 ^ 2 / t70 ^ 2 * t69);
	t50 = 0.1e1 / t53;
	t49 = 0.1e1 / t77;
	t46 = 0.1e1 / t48;
	t45 = 0.1e1 / (0.1e1 + t83);
	t44 = (t58 * t69 * t71 + t60 * t68) * t67 * t57;
	t43 = t77 * t49;
	t1 = [t61 * t57 * t82, t44, 0, 0, 0, 0; (-t58 * t46 - (-t54 + (-t82 * t84 + t54) * t57) * t83) * t45, (t62 * t46 - (t55 * t70 * t71 - t54 * t60 + (t54 * t80 - t84) * t44) * t61 * t47) * t45, 0, 0, 0, 0; ((-t60 * t64 - t65 * t79) * t50 - (-t60 * t65 + t64 * t79) * t85) * t49, (-t64 * t50 + t65 * t85) * t61 * t49, 0, t43, t43, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (1756->38), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
	t101 = pkin(12) + qJ(4) + qJ(5);
	t100 = cos(t101);
	t102 = sin(pkin(6));
	t109 = cos(qJ(1));
	t116 = t102 * t109;
	t103 = cos(pkin(6));
	t105 = sin(qJ(2));
	t113 = t109 * t105;
	t106 = sin(qJ(1));
	t108 = cos(qJ(2));
	t114 = t106 * t108;
	t94 = t103 * t113 + t114;
	t99 = sin(t101);
	t83 = t100 * t116 + t94 * t99;
	t119 = t102 * t105;
	t91 = -t103 * t100 + t99 * t119;
	t78 = atan2(-t83, t91);
	t75 = sin(t78);
	t76 = cos(t78);
	t73 = -t75 * t83 + t76 * t91;
	t72 = 0.1e1 / t73 ^ 2;
	t118 = t102 * t106;
	t112 = t109 * t108;
	t115 = t106 * t105;
	t96 = -t103 * t115 + t112;
	t87 = -t100 * t118 + t96 * t99;
	t126 = t72 * t87;
	t125 = t76 * t83;
	t107 = cos(qJ(6));
	t104 = sin(qJ(6));
	t95 = t103 * t114 + t113;
	t121 = t95 * t104;
	t88 = t96 * t100 + t99 * t118;
	t82 = t88 * t107 + t121;
	t80 = 0.1e1 / t82 ^ 2;
	t120 = t95 * t107;
	t81 = t88 * t104 - t120;
	t124 = t80 * t81;
	t90 = 0.1e1 / t91 ^ 2;
	t123 = t83 * t90;
	t122 = t87 ^ 2 * t72;
	t117 = t102 * t108;
	t111 = t81 ^ 2 * t80 + 0.1e1;
	t85 = t94 * t100 - t99 * t116;
	t110 = -t75 * t91 - t125;
	t93 = t103 * t112 - t115;
	t92 = t100 * t119 + t103 * t99;
	t89 = 0.1e1 / t91;
	t79 = 0.1e1 / t82;
	t77 = 0.1e1 / (t83 ^ 2 * t90 + 0.1e1);
	t74 = 0.1e1 / t111;
	t71 = 0.1e1 / t73;
	t70 = 0.1e1 / (0.1e1 + t122);
	t69 = (t117 * t123 - t89 * t93) * t99 * t77;
	t68 = (t92 * t123 - t85 * t89) * t77;
	t67 = (-t104 * t79 + t107 * t124) * t87 * t74;
	t66 = (t88 * t71 - (t110 * t68 - t75 * t85 + t76 * t92) * t126) * t70;
	t1 = [-t87 * t89 * t77, t69, 0, t68, t68, 0; (-t83 * t71 - (-t75 + (t89 * t125 + t75) * t77) * t122) * t70, (-t95 * t99 * t71 - ((t76 * t117 - t75 * t93) * t99 + t110 * t69) * t126) * t70, 0, t66, t66, 0; ((-t104 * t85 - t93 * t107) * t79 - (t93 * t104 - t107 * t85) * t124) * t74, ((-t100 * t121 - t96 * t107) * t79 - (-t100 * t120 + t96 * t104) * t124) * t74, 0, t67, t67, t111 * t74;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end