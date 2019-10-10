% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR4
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
%   Wie in S6PRRPRR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRPRR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
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
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (256->25), mult. (739->65), div. (57->9), fcn. (1061->11), ass. (0->42)
	t55 = sin(pkin(11));
	t57 = cos(pkin(11));
	t62 = cos(qJ(2));
	t58 = cos(pkin(6));
	t60 = sin(qJ(2));
	t65 = t58 * t60;
	t50 = t55 * t62 + t57 * t65;
	t59 = sin(qJ(3));
	t56 = sin(pkin(6));
	t61 = cos(qJ(3));
	t67 = t56 * t61;
	t41 = t50 * t59 + t57 * t67;
	t68 = t56 * t59;
	t53 = -t58 * t61 + t60 * t68;
	t39 = atan2(-t41, t53);
	t35 = sin(t39);
	t36 = cos(t39);
	t34 = -t35 * t41 + t36 * t53;
	t33 = 0.1e1 / t34 ^ 2;
	t52 = -t55 * t65 + t57 * t62;
	t44 = t52 * t59 - t55 * t67;
	t71 = t33 * t44;
	t48 = 0.1e1 / t53 ^ 2;
	t70 = t41 * t48;
	t45 = t52 * t61 + t55 * t68;
	t40 = 0.1e1 / t45 ^ 2;
	t64 = t58 * t62;
	t51 = -t55 * t64 - t57 * t60;
	t69 = t51 ^ 2 * t40;
	t66 = t56 * t62;
	t63 = -t35 * t53 - t36 * t41;
	t54 = t58 * t59 + t60 * t67;
	t49 = -t55 * t60 + t57 * t64;
	t47 = 0.1e1 / t53;
	t43 = t50 * t61 - t57 * t68;
	t38 = 0.1e1 / (t41 ^ 2 * t48 + 0.1e1);
	t37 = 0.1e1 / (0.1e1 + t69);
	t32 = 0.1e1 / t34;
	t31 = 0.1e1 / (t44 ^ 2 * t33 + 0.1e1);
	t30 = (-t47 * t49 + t66 * t70) * t59 * t38;
	t29 = (-t43 * t47 + t54 * t70) * t38;
	t1 = [0, t30, t29, 0, 0, 0; 0, (t51 * t59 * t32 - ((-t35 * t49 + t36 * t66) * t59 + t63 * t30) * t71) * t31, (t45 * t32 - (t63 * t29 - t35 * t43 + t36 * t54) * t71) * t31, 0, 0, 0; 0, (-t52 / t45 - t61 * t69) * t37, t44 * t51 * t40 * t37, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (183->21), mult. (524->53), div. (53->11), fcn. (745->13), ass. (0->35)
	t55 = sin(pkin(11));
	t57 = cos(pkin(11));
	t64 = cos(qJ(2));
	t58 = cos(pkin(6));
	t61 = sin(qJ(2));
	t67 = t58 * t61;
	t50 = -t55 * t67 + t57 * t64;
	t60 = sin(qJ(3));
	t63 = cos(qJ(3));
	t56 = sin(pkin(6));
	t69 = t55 * t56;
	t41 = t50 * t60 - t63 * t69;
	t42 = t50 * t63 + t60 * t69;
	t59 = sin(qJ(5));
	t62 = cos(qJ(5));
	t40 = t41 * t59 + t42 * t62;
	t38 = 0.1e1 / t40 ^ 2;
	t39 = -t41 * t62 + t42 * t59;
	t70 = t39 ^ 2 * t38;
	t68 = t56 * t64;
	t66 = t58 * t64;
	t65 = 0.1e1 + t70;
	t54 = 0.1e1 / t64 ^ 2;
	t49 = -t55 * t66 - t57 * t61;
	t48 = t55 * t64 + t57 * t67;
	t47 = t55 * t61 - t57 * t66;
	t46 = atan2(t47, t68);
	t44 = cos(t46);
	t43 = sin(t46);
	t37 = 0.1e1 / t40;
	t36 = t43 * t47 + t44 * t68;
	t35 = 0.1e1 / t36 ^ 2;
	t33 = (t48 / t64 + t61 * t47 * t54) / t56 / (0.1e1 + t47 ^ 2 / t56 ^ 2 * t54);
	t32 = 0.1e1 / t65;
	t1 = [0, t33, 0, 0, 0, 0; 0, (-t50 / t36 - (-t44 * t56 * t61 + t43 * t48 + (-t43 * t68 + t44 * t47) * t33) * t49 * t35) / (t49 ^ 2 * t35 + 0.1e1), 0, 0, 0, 0; 0, ((t59 * t63 - t60 * t62) * t37 - (t59 * t60 + t62 * t63) * t39 * t38) * t32 * t49, (-t40 * t37 - t70) * t32, 0, t65 * t32, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:30
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (1084->41), mult. (2952->94), div. (95->9), fcn. (4136->15), ass. (0->60)
	t113 = cos(pkin(11));
	t121 = cos(qJ(3));
	t111 = sin(pkin(11));
	t122 = cos(qJ(2));
	t114 = cos(pkin(6));
	t118 = sin(qJ(2));
	t130 = t114 * t118;
	t124 = t111 * t122 + t113 * t130;
	t112 = sin(pkin(6));
	t117 = sin(qJ(3));
	t132 = t112 * t117;
	t100 = -t113 * t132 + t121 * t124;
	t116 = sin(qJ(5));
	t120 = cos(qJ(5));
	t131 = t112 * t121;
	t123 = t113 * t131 + t117 * t124;
	t85 = t100 * t116 - t120 * t123;
	t107 = -t114 * t121 + t118 * t132;
	t108 = t114 * t117 + t118 * t131;
	t96 = -t107 * t120 + t108 * t116;
	t79 = atan2(-t85, t96);
	t76 = sin(t79);
	t77 = cos(t79);
	t127 = -t76 * t96 - t77 * t85;
	t74 = -t76 * t85 + t77 * t96;
	t73 = 0.1e1 / t74 ^ 2;
	t106 = -t111 * t130 + t113 * t122;
	t101 = t106 * t117 - t111 * t131;
	t102 = t106 * t121 + t111 * t132;
	t89 = -t101 * t120 + t102 * t116;
	t137 = t73 * t89;
	t94 = 0.1e1 / t96 ^ 2;
	t135 = t85 * t94;
	t78 = 0.1e1 / (t85 ^ 2 * t94 + 0.1e1);
	t87 = t100 * t120 + t116 * t123;
	t93 = 0.1e1 / t96;
	t97 = t107 * t116 + t108 * t120;
	t143 = (-t97 * t135 + t87 * t93) * t78;
	t71 = 0.1e1 / (t73 * t89 ^ 2 + 0.1e1);
	t72 = 0.1e1 / t74;
	t90 = t101 * t116 + t102 * t120;
	t146 = ((t127 * t143 + t76 * t87 - t77 * t97) * t137 + t90 * t72) * t71;
	t115 = sin(qJ(6));
	t119 = cos(qJ(6));
	t129 = t114 * t122;
	t105 = -t111 * t129 - t113 * t118;
	t83 = t105 * t115 + t119 * t90;
	t81 = 0.1e1 / t83 ^ 2;
	t82 = -t105 * t119 + t115 * t90;
	t136 = t81 * t82;
	t128 = t81 * t82 ^ 2 + 0.1e1;
	t75 = 0.1e1 / t128;
	t80 = 0.1e1 / t83;
	t138 = (-t115 * t80 + t119 * t136) * t89 * t75;
	t126 = t116 * t121 - t117 * t120;
	t103 = t126 * t122 * t112;
	t92 = (t116 * t117 + t120 * t121) * t105;
	t91 = t126 * (-t111 * t118 + t113 * t129);
	t70 = (t103 * t135 - t91 * t93) * t78;
	t1 = [0, t70, t143, 0, -t143, 0; 0, (-(t77 * t103 + t127 * t70 - t76 * t91) * t137 + t126 * t72 * t105) * t71, -t146, 0, t146, 0; 0, ((t106 * t119 + t115 * t92) * t80 - (-t106 * t115 + t119 * t92) * t136) * t75, -t138, 0, t138, t128 * t75;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end