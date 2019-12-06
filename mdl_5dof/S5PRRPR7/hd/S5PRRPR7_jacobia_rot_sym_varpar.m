% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5PRRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRRPR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:32
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (101->18), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->29)
	t41 = sin(pkin(9));
	t42 = sin(pkin(5));
	t53 = t41 * t42;
	t48 = cos(qJ(2));
	t52 = t42 * t48;
	t44 = cos(pkin(5));
	t46 = sin(qJ(2));
	t51 = t44 * t46;
	t50 = t44 * t48;
	t43 = cos(pkin(9));
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
	t1 = [0, t21, 0, 0, 0; 0, (t37 / t24 - (t30 * t42 * t46 - t29 * t35 + (t29 * t52 - t30 * t33) * t21) * t36 * t23) / (t36 ^ 2 * t23 + 0.1e1), 0, 0, 0; 0, (-t45 / t28 + t47 * t27 * t26) * t36 * t25, t49 * t25, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:32
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (303->30), mult. (866->78), div. (60->9), fcn. (1237->13), ass. (0->48)
	t61 = sin(pkin(9));
	t64 = cos(pkin(9));
	t69 = cos(qJ(2));
	t65 = cos(pkin(5));
	t67 = sin(qJ(2));
	t72 = t65 * t67;
	t54 = t61 * t69 + t64 * t72;
	t66 = sin(qJ(3));
	t62 = sin(pkin(5));
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
	t60 = sin(pkin(10));
	t63 = cos(pkin(10));
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
	t1 = [0, t32, t31, 0, 0; 0, (-t55 * t66 * t34 - ((-t42 * t53 + t43 * t73) * t66 + t70 * t32) * t79) * t33, (t50 * t34 - (t70 * t31 - t42 * t48 + t43 * t58) * t79) * t33, 0, 0; 0, ((-t56 * t63 - t60 * t76) * t38 - (t56 * t60 - t63 * t76) * t78) * t37, (-t38 * t60 + t63 * t78) * t49 * t37, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:32
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (595->40), mult. (1655->100), div. (65->9), fcn. (2316->15), ass. (0->59)
	t86 = sin(pkin(5));
	t91 = sin(qJ(3));
	t102 = t86 * t91;
	t89 = cos(pkin(5));
	t92 = sin(qJ(2));
	t100 = t89 * t92;
	t85 = sin(pkin(9));
	t88 = cos(pkin(9));
	t95 = cos(qJ(2));
	t79 = t88 * t100 + t85 * t95;
	t84 = sin(pkin(10));
	t87 = cos(pkin(10));
	t94 = cos(qJ(3));
	t99 = t89 * t95;
	t96 = -t85 * t92 + t88 * t99;
	t65 = (-t88 * t102 + t79 * t94) * t84 + t96 * t87;
	t101 = t86 * t94;
	t73 = (t92 * t101 + t89 * t91) * t84 + t86 * t95 * t87;
	t64 = atan2(-t65, t73);
	t61 = sin(t64);
	t62 = cos(t64);
	t55 = -t61 * t65 + t62 * t73;
	t54 = 0.1e1 / t55 ^ 2;
	t80 = t85 * t99 + t88 * t92;
	t105 = t80 * t87;
	t81 = -t85 * t100 + t88 * t95;
	t76 = t85 * t102 + t81 * t94;
	t67 = t76 * t84 - t105;
	t110 = t54 * t67;
	t75 = -t85 * t101 + t81 * t91;
	t90 = sin(qJ(5));
	t107 = t75 * t90;
	t68 = t76 * t87 + t80 * t84;
	t93 = cos(qJ(5));
	t60 = t68 * t93 + t107;
	t58 = 0.1e1 / t60 ^ 2;
	t106 = t75 * t93;
	t59 = t68 * t90 - t106;
	t109 = t58 * t59;
	t72 = 0.1e1 / t73 ^ 2;
	t108 = t65 * t72;
	t104 = t80 * t91;
	t103 = t84 * t94;
	t98 = t58 * t59 ^ 2 + 0.1e1;
	t97 = -t61 * t73 - t62 * t65;
	t82 = -t92 * t102 + t89 * t94;
	t77 = (t95 * t103 - t87 * t92) * t86;
	t74 = -t88 * t101 - t79 * t91;
	t71 = 0.1e1 / t73;
	t70 = -t94 * t105 + t81 * t84;
	t69 = t96 * t103 - t79 * t87;
	t63 = 0.1e1 / (t65 ^ 2 * t72 + 0.1e1);
	t57 = 0.1e1 / t60;
	t56 = 0.1e1 / t98;
	t53 = 0.1e1 / t55;
	t52 = 0.1e1 / (t54 * t67 ^ 2 + 0.1e1);
	t51 = (t82 * t108 - t71 * t74) * t84 * t63;
	t50 = (t77 * t108 - t69 * t71) * t63;
	t1 = [0, t50, t51, 0, 0; 0, ((-t80 * t103 - t81 * t87) * t53 - (t97 * t50 - t61 * t69 + t62 * t77) * t110) * t52, (-t75 * t84 * t53 - ((-t61 * t74 + t62 * t82) * t84 + t97 * t51) * t110) * t52, 0, 0; 0, ((t93 * t104 + t70 * t90) * t57 - (-t90 * t104 + t70 * t93) * t109) * t56, ((-t87 * t107 - t76 * t93) * t57 - (-t87 * t106 + t76 * t90) * t109) * t56, 0, t98 * t56;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end