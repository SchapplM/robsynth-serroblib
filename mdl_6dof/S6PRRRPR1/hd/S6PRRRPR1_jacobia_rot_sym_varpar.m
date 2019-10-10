% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR1
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
%   Wie in S6PRRRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRPR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:43
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
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (206->19), mult. (316->43), div. (52->11), fcn. (458->11), ass. (0->31)
	t64 = sin(pkin(11));
	t65 = sin(pkin(6));
	t74 = t64 * t65;
	t69 = cos(qJ(2));
	t73 = t65 * t69;
	t67 = cos(pkin(6));
	t68 = sin(qJ(2));
	t72 = t67 * t68;
	t71 = t67 * t69;
	t66 = cos(pkin(11));
	t57 = -t64 * t72 + t66 * t69;
	t61 = qJ(3) + qJ(4) + pkin(12);
	t59 = sin(t61);
	t60 = cos(t61);
	t48 = t57 * t60 + t59 * t74;
	t46 = 0.1e1 / t48 ^ 2;
	t47 = t57 * t59 - t60 * t74;
	t70 = t47 ^ 2 * t46 + 0.1e1;
	t63 = 0.1e1 / t69 ^ 2;
	t56 = t64 * t71 + t66 * t68;
	t55 = t64 * t69 + t66 * t72;
	t53 = t64 * t68 - t66 * t71;
	t51 = atan2(-t53, -t73);
	t50 = cos(t51);
	t49 = sin(t51);
	t45 = 0.1e1 / t70;
	t44 = -t49 * t53 - t50 * t73;
	t43 = 0.1e1 / t44 ^ 2;
	t41 = (t55 / t69 + t68 * t53 * t63) / t65 / (0.1e1 + t53 ^ 2 / t65 ^ 2 * t63);
	t40 = t70 * t45;
	t1 = [0, t41, 0, 0, 0, 0; 0, (t57 / t44 - (t50 * t65 * t68 - t49 * t55 + (t49 * t73 - t50 * t53) * t41) * t56 * t43) / (t56 ^ 2 * t43 + 0.1e1), 0, 0, 0, 0; 0, (-t59 / t48 + t60 * t47 * t46) * t56 * t45, t40, t40, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:43
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (1407->31), mult. (1387->78), div. (95->9), fcn. (1976->13), ass. (0->55)
	t95 = sin(pkin(6));
	t96 = cos(pkin(11));
	t110 = t95 * t96;
	t101 = cos(qJ(2));
	t94 = sin(pkin(11));
	t105 = t94 * t101;
	t97 = cos(pkin(6));
	t99 = sin(qJ(2));
	t108 = t97 * t99;
	t87 = t96 * t108 + t105;
	t93 = qJ(3) + qJ(4) + pkin(12);
	t91 = sin(t93);
	t92 = cos(t93);
	t77 = t92 * t110 + t87 * t91;
	t109 = t95 * t99;
	t84 = t91 * t109 - t97 * t92;
	t72 = atan2(-t77, t84);
	t69 = sin(t72);
	t70 = cos(t72);
	t67 = -t69 * t77 + t70 * t84;
	t66 = 0.1e1 / t67 ^ 2;
	t111 = t94 * t95;
	t104 = t96 * t101;
	t89 = -t94 * t108 + t104;
	t80 = -t92 * t111 + t89 * t91;
	t115 = t66 * t80;
	t100 = cos(qJ(6));
	t88 = t97 * t105 + t96 * t99;
	t98 = sin(qJ(6));
	t112 = t88 * t98;
	t81 = t91 * t111 + t89 * t92;
	t76 = t81 * t100 + t112;
	t74 = 0.1e1 / t76 ^ 2;
	t106 = t88 * t100;
	t75 = t81 * t98 - t106;
	t114 = t74 * t75;
	t83 = 0.1e1 / t84 ^ 2;
	t113 = t77 * t83;
	t107 = t101 * t95;
	t103 = t75 ^ 2 * t74 + 0.1e1;
	t102 = -t69 * t84 - t70 * t77;
	t86 = t97 * t104 - t94 * t99;
	t85 = t92 * t109 + t97 * t91;
	t82 = 0.1e1 / t84;
	t79 = -t91 * t110 + t87 * t92;
	t73 = 0.1e1 / t76;
	t71 = 0.1e1 / (t77 ^ 2 * t83 + 0.1e1);
	t68 = 0.1e1 / t103;
	t65 = 0.1e1 / t67;
	t64 = 0.1e1 / (t80 ^ 2 * t66 + 0.1e1);
	t63 = (t107 * t113 - t82 * t86) * t91 * t71;
	t62 = (t85 * t113 - t79 * t82) * t71;
	t61 = (t100 * t114 - t73 * t98) * t80 * t68;
	t60 = (t81 * t65 - (t102 * t62 - t69 * t79 + t70 * t85) * t115) * t64;
	t1 = [0, t63, t62, t62, 0, 0; 0, (-t88 * t91 * t65 - ((t70 * t107 - t69 * t86) * t91 + t102 * t63) * t115) * t64, t60, t60, 0, 0; 0, ((-t89 * t100 - t92 * t112) * t73 - (-t92 * t106 + t89 * t98) * t114) * t68, t61, t61, 0, t103 * t68;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end