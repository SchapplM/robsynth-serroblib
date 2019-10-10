% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR1
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
%   Wie in S6PRRRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (101->18), mult. (268->43), div. (47->11), fcn. (395->11), ass. (0->29)
	t41 = sin(pkin(12));
	t42 = sin(pkin(6));
	t53 = t41 * t42;
	t48 = cos(qJ(2));
	t52 = t42 * t48;
	t44 = cos(pkin(6));
	t46 = sin(qJ(2));
	t51 = t44 * t46;
	t50 = t44 * t48;
	t43 = cos(pkin(12));
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
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (162->19), mult. (316->43), div. (52->11), fcn. (458->11), ass. (0->31)
	t59 = sin(pkin(12));
	t60 = sin(pkin(6));
	t69 = t59 * t60;
	t64 = cos(qJ(2));
	t68 = t60 * t64;
	t62 = cos(pkin(6));
	t63 = sin(qJ(2));
	t67 = t62 * t63;
	t66 = t62 * t64;
	t61 = cos(pkin(12));
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
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (251->19), mult. (364->43), div. (57->11), fcn. (521->11), ass. (0->31)
	t63 = sin(pkin(12));
	t64 = sin(pkin(6));
	t73 = t63 * t64;
	t68 = cos(qJ(2));
	t72 = t64 * t68;
	t66 = cos(pkin(6));
	t67 = sin(qJ(2));
	t71 = t66 * t67;
	t70 = t66 * t68;
	t65 = cos(pkin(12));
	t56 = -t63 * t71 + t65 * t68;
	t60 = qJ(3) + qJ(4) + qJ(5);
	t58 = sin(t60);
	t59 = cos(t60);
	t47 = t56 * t59 + t58 * t73;
	t45 = 0.1e1 / t47 ^ 2;
	t46 = t56 * t58 - t59 * t73;
	t69 = t46 ^ 2 * t45 + 0.1e1;
	t62 = 0.1e1 / t68 ^ 2;
	t55 = t63 * t70 + t65 * t67;
	t54 = t63 * t68 + t65 * t71;
	t52 = t63 * t67 - t65 * t70;
	t50 = atan2(-t52, -t72);
	t49 = cos(t50);
	t48 = sin(t50);
	t44 = 0.1e1 / t69;
	t43 = -t48 * t52 - t49 * t72;
	t42 = 0.1e1 / t43 ^ 2;
	t40 = (t54 / t68 + t67 * t52 * t62) / t64 / (0.1e1 + t52 ^ 2 / t64 ^ 2 * t62);
	t39 = t69 * t44;
	t1 = [0, t40, 0, 0, 0, 0; 0, (t56 / t43 - (t49 * t64 * t67 - t48 * t54 + (t48 * t72 - t49 * t52) * t40) * t55 * t42) / (t55 ^ 2 * t42 + 0.1e1), 0, 0, 0, 0; 0, (-t58 / t47 + t59 * t46 * t45) * t55 * t44, t39, t39, t39, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (1866->31), mult. (1825->77), div. (125->9), fcn. (2603->13), ass. (0->56)
	t97 = sin(pkin(6));
	t98 = cos(pkin(12));
	t114 = t97 * t98;
	t101 = sin(qJ(2));
	t107 = t98 * t101;
	t103 = cos(qJ(2));
	t96 = sin(pkin(12));
	t108 = t96 * t103;
	t99 = cos(pkin(6));
	t89 = t99 * t107 + t108;
	t95 = qJ(3) + qJ(4) + qJ(5);
	t93 = sin(t95);
	t94 = cos(t95);
	t79 = t94 * t114 + t89 * t93;
	t113 = t101 * t97;
	t86 = t93 * t113 - t99 * t94;
	t76 = atan2(-t79, t86);
	t71 = sin(t76);
	t72 = cos(t76);
	t69 = -t71 * t79 + t72 * t86;
	t68 = 0.1e1 / t69 ^ 2;
	t115 = t96 * t97;
	t106 = t98 * t103;
	t109 = t96 * t101;
	t91 = -t99 * t109 + t106;
	t82 = -t94 * t115 + t91 * t93;
	t118 = t68 * t82;
	t102 = cos(qJ(6));
	t100 = sin(qJ(6));
	t90 = t99 * t108 + t107;
	t111 = t90 * t100;
	t83 = t93 * t115 + t91 * t94;
	t78 = t83 * t102 + t111;
	t75 = 0.1e1 / t78 ^ 2;
	t110 = t90 * t102;
	t77 = t83 * t100 - t110;
	t117 = t75 * t77;
	t85 = 0.1e1 / t86 ^ 2;
	t116 = t79 * t85;
	t112 = t103 * t97;
	t105 = t77 ^ 2 * t75 + 0.1e1;
	t104 = -t71 * t86 - t72 * t79;
	t88 = t99 * t106 - t109;
	t87 = t94 * t113 + t99 * t93;
	t84 = 0.1e1 / t86;
	t81 = -t93 * t114 + t89 * t94;
	t74 = 0.1e1 / t78;
	t73 = 0.1e1 / (t79 ^ 2 * t85 + 0.1e1);
	t70 = 0.1e1 / t105;
	t67 = 0.1e1 / t69;
	t66 = 0.1e1 / (t82 ^ 2 * t68 + 0.1e1);
	t65 = (t112 * t116 - t84 * t88) * t93 * t73;
	t64 = (t87 * t116 - t81 * t84) * t73;
	t63 = (-t100 * t74 + t102 * t117) * t82 * t70;
	t62 = (t83 * t67 - (t104 * t64 - t71 * t81 + t72 * t87) * t118) * t66;
	t1 = [0, t65, t64, t64, t64, 0; 0, (-t90 * t93 * t67 - ((t72 * t112 - t71 * t88) * t93 + t104 * t65) * t118) * t66, t62, t62, t62, 0; 0, ((-t91 * t102 - t94 * t111) * t74 - (t91 * t100 - t94 * t110) * t117) * t70, t63, t63, t63, t105 * t70;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end