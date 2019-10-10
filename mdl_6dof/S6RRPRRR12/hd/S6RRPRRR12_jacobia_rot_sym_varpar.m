% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR12
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
%   Wie in S6RRPRRR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRR12_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR12_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:01
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:02
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
	% StartTime: 2019-10-10 11:07:01
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (126->20), mult. (316->53), div. (70->11), fcn. (497->9), ass. (0->33)
	t43 = sin(qJ(2));
	t44 = sin(qJ(1));
	t45 = cos(qJ(2));
	t46 = cos(qJ(1));
	t49 = cos(pkin(6));
	t47 = t46 * t49;
	t30 = t44 * t43 - t45 * t47;
	t42 = sin(pkin(6));
	t50 = t42 * t45;
	t27 = atan2(-t30, -t50);
	t26 = cos(t27);
	t54 = t26 * t30;
	t48 = t44 * t49;
	t34 = -t43 * t48 + t46 * t45;
	t36 = 0.1e1 / t42;
	t37 = 0.1e1 / t42 ^ 2;
	t39 = 0.1e1 / t44 ^ 2;
	t53 = 0.1e1 / (t34 ^ 2 * t39 * t37 + 0.1e1) * t36;
	t25 = sin(t27);
	t24 = -t25 * t30 - t26 * t50;
	t23 = 0.1e1 / t24 ^ 2;
	t33 = t46 * t43 + t45 * t48;
	t52 = t33 ^ 2 * t23;
	t40 = 0.1e1 / t45;
	t51 = t36 * t40;
	t41 = 0.1e1 / t45 ^ 2;
	t38 = 0.1e1 / t44;
	t32 = t43 * t47 + t44 * t45;
	t28 = 0.1e1 / (t30 ^ 2 * t37 * t41 + 0.1e1);
	t22 = 0.1e1 / t24;
	t21 = 0.1e1 / (0.1e1 + t52);
	t20 = (t30 * t41 * t43 + t32 * t40) * t36 * t28;
	t1 = [t33 * t28 * t51, t20, 0, 0, 0, 0; (-t30 * t22 - (-t25 + (-t51 * t54 + t25) * t28) * t52) * t21, (t34 * t22 - (t26 * t42 * t43 - t25 * t32 + (t25 * t50 - t54) * t20) * t33 * t23) * t21, 0, 0, 0, 0; (-t34 * t39 * t46 - t32 * t38) * t53, -t33 * t38 * t53, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (149->24), mult. (451->64), div. (72->11), fcn. (673->11), ass. (0->40)
	t54 = sin(qJ(2));
	t57 = cos(qJ(2));
	t58 = cos(qJ(1));
	t55 = sin(qJ(1));
	t62 = cos(pkin(6));
	t60 = t55 * t62;
	t45 = t58 * t54 + t57 * t60;
	t53 = sin(qJ(4));
	t56 = cos(qJ(4));
	t52 = sin(pkin(6));
	t64 = t52 * t55;
	t37 = t45 * t53 + t56 * t64;
	t35 = 0.1e1 / t37 ^ 2;
	t36 = -t45 * t56 + t53 * t64;
	t69 = t35 * t36;
	t59 = t58 * t62;
	t43 = t54 * t59 + t55 * t57;
	t65 = t52 * t54;
	t41 = atan2(-t43, t65);
	t39 = cos(t41);
	t68 = t39 * t43;
	t38 = sin(t41);
	t32 = -t38 * t43 + t39 * t65;
	t31 = 0.1e1 / t32 ^ 2;
	t46 = -t54 * t60 + t58 * t57;
	t67 = t46 ^ 2 * t31;
	t49 = 0.1e1 / t52;
	t50 = 0.1e1 / t54;
	t66 = t49 * t50;
	t63 = t52 * t58;
	t61 = t36 ^ 2 * t35 + 0.1e1;
	t51 = 0.1e1 / t54 ^ 2;
	t42 = t55 * t54 - t57 * t59;
	t40 = 0.1e1 / (0.1e1 + t43 ^ 2 / t52 ^ 2 * t51);
	t34 = 0.1e1 / t37;
	t33 = 0.1e1 / t61;
	t30 = 0.1e1 / t32;
	t29 = 0.1e1 / (0.1e1 + t67);
	t28 = (t43 * t51 * t57 + t42 * t50) * t49 * t40;
	t1 = [-t46 * t40 * t66, t28, 0, 0, 0, 0; (-t43 * t30 - (-t38 + (t66 * t68 + t38) * t40) * t67) * t29, (-t45 * t30 - (t39 * t52 * t57 + t38 * t42 + (-t38 * t65 - t68) * t28) * t46 * t31) * t29, 0, 0, 0, 0; ((t42 * t56 + t53 * t63) * t34 - (-t42 * t53 + t56 * t63) * t69) * t33, (-t56 * t34 - t53 * t69) * t46 * t33, 0, t61 * t33, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (228->25), mult. (499->64), div. (77->11), fcn. (736->11), ass. (0->42)
	t72 = sin(qJ(2));
	t74 = cos(qJ(2));
	t75 = cos(qJ(1));
	t73 = sin(qJ(1));
	t79 = cos(pkin(6));
	t77 = t73 * t79;
	t61 = t75 * t72 + t74 * t77;
	t70 = qJ(4) + qJ(5);
	t65 = sin(t70);
	t66 = cos(t70);
	t71 = sin(pkin(6));
	t81 = t71 * t73;
	t53 = t61 * t65 + t66 * t81;
	t51 = 0.1e1 / t53 ^ 2;
	t52 = -t61 * t66 + t65 * t81;
	t86 = t51 * t52;
	t76 = t75 * t79;
	t59 = t72 * t76 + t73 * t74;
	t82 = t71 * t72;
	t57 = atan2(-t59, t82);
	t55 = cos(t57);
	t85 = t55 * t59;
	t54 = sin(t57);
	t48 = -t54 * t59 + t55 * t82;
	t47 = 0.1e1 / t48 ^ 2;
	t62 = -t72 * t77 + t75 * t74;
	t84 = t62 ^ 2 * t47;
	t67 = 0.1e1 / t71;
	t68 = 0.1e1 / t72;
	t83 = t67 * t68;
	t80 = t71 * t75;
	t78 = t52 ^ 2 * t51 + 0.1e1;
	t69 = 0.1e1 / t72 ^ 2;
	t58 = t73 * t72 - t74 * t76;
	t56 = 0.1e1 / (0.1e1 + t59 ^ 2 / t71 ^ 2 * t69);
	t50 = 0.1e1 / t53;
	t49 = 0.1e1 / t78;
	t46 = 0.1e1 / t48;
	t45 = 0.1e1 / (0.1e1 + t84);
	t44 = (t59 * t69 * t74 + t58 * t68) * t67 * t56;
	t43 = t78 * t49;
	t1 = [-t62 * t56 * t83, t44, 0, 0, 0, 0; (-t59 * t46 - (-t54 + (t83 * t85 + t54) * t56) * t84) * t45, (-t61 * t46 - (t55 * t71 * t74 + t54 * t58 + (-t54 * t82 - t85) * t44) * t62 * t47) * t45, 0, 0, 0, 0; ((t58 * t66 + t65 * t80) * t50 - (-t58 * t65 + t66 * t80) * t86) * t49, (-t66 * t50 - t65 * t86) * t62 * t49, 0, t43, t43, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:07:02
	% EndTime: 2019-10-10 11:07:02
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (1185->37), mult. (1738->90), div. (115->9), fcn. (2485->13), ass. (0->58)
	t103 = cos(qJ(1));
	t96 = sin(pkin(6));
	t112 = t103 * t96;
	t102 = cos(qJ(2));
	t108 = t103 * t102;
	t100 = sin(qJ(1));
	t99 = sin(qJ(2));
	t114 = t100 * t99;
	t97 = cos(pkin(6));
	t89 = -t97 * t108 + t114;
	t95 = qJ(4) + qJ(5);
	t93 = sin(t95);
	t94 = cos(t95);
	t81 = t93 * t112 + t89 * t94;
	t113 = t102 * t96;
	t87 = t94 * t113 + t97 * t93;
	t78 = atan2(t81, t87);
	t73 = sin(t78);
	t74 = cos(t78);
	t69 = t73 * t81 + t74 * t87;
	t68 = 0.1e1 / t69 ^ 2;
	t109 = t100 * t102;
	t111 = t103 * t99;
	t104 = t97 * t109 + t111;
	t115 = t100 * t96;
	t79 = -t104 * t94 + t93 * t115;
	t122 = t68 * t79;
	t101 = cos(qJ(6));
	t91 = -t97 * t114 + t108;
	t98 = sin(qJ(6));
	t117 = t91 * t98;
	t80 = t104 * t93 + t94 * t115;
	t76 = t80 * t101 + t117;
	t72 = 0.1e1 / t76 ^ 2;
	t110 = t91 * t101;
	t75 = t80 * t98 - t110;
	t121 = t72 * t75;
	t120 = t74 * t81;
	t119 = t79 ^ 2 * t68;
	t86 = 0.1e1 / t87 ^ 2;
	t118 = t81 * t86;
	t116 = t96 * t99;
	t107 = t75 ^ 2 * t72 + 0.1e1;
	t106 = -t73 * t87 + t120;
	t105 = t94 * t112 - t89 * t93;
	t90 = t97 * t111 + t109;
	t88 = -t93 * t113 + t97 * t94;
	t85 = 0.1e1 / t87;
	t77 = 0.1e1 / (t81 ^ 2 * t86 + 0.1e1);
	t71 = 0.1e1 / t76;
	t70 = 0.1e1 / t107;
	t67 = 0.1e1 / t69;
	t66 = 0.1e1 / (0.1e1 + t119);
	t65 = (t116 * t118 + t85 * t90) * t94 * t77;
	t64 = (t105 * t85 - t88 * t118) * t77;
	t63 = (t101 * t121 - t98 * t71) * t79 * t70;
	t62 = (t80 * t67 - (t105 * t73 + t106 * t64 + t74 * t88) * t122) * t66;
	t1 = [-t79 * t85 * t77, t65, 0, t64, t64, 0; (t81 * t67 - (-t73 + (-t85 * t120 + t73) * t77) * t119) * t66, (-t91 * t94 * t67 - ((-t74 * t116 + t73 * t90) * t94 + t106 * t65) * t122) * t66, 0, t62, t62, 0; ((t90 * t101 + t105 * t98) * t71 - (t101 * t105 - t90 * t98) * t121) * t70, ((t104 * t101 + t93 * t117) * t71 - (-t104 * t98 + t93 * t110) * t121) * t70, 0, t63, t63, t107 * t70;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end