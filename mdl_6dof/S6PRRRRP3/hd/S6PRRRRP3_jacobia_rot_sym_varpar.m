% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRP3
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
%   Wie in S6PRRRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:05
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRRRRP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:05:49
	% EndTime: 2019-10-09 23:05:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:05:49
	% EndTime: 2019-10-09 23:05:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:05:49
	% EndTime: 2019-10-09 23:05:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:05:49
	% EndTime: 2019-10-09 23:05:49
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
	% StartTime: 2019-10-09 23:05:49
	% EndTime: 2019-10-09 23:05:49
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (334->30), mult. (949->79), div. (65->9), fcn. (1349->13), ass. (0->49)
	t62 = sin(pkin(11));
	t64 = cos(pkin(11));
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
	% StartTime: 2019-10-09 23:05:49
	% EndTime: 2019-10-09 23:05:49
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (427->31), mult. (1032->79), div. (70->9), fcn. (1461->13), ass. (0->51)
	t80 = sin(pkin(11));
	t82 = cos(pkin(11));
	t87 = cos(qJ(2));
	t83 = cos(pkin(6));
	t85 = sin(qJ(2));
	t91 = t83 * t85;
	t71 = t80 * t87 + t82 * t91;
	t84 = sin(qJ(3));
	t81 = sin(pkin(6));
	t86 = cos(qJ(3));
	t93 = t81 * t86;
	t63 = t71 * t84 + t82 * t93;
	t94 = t81 * t84;
	t74 = -t83 * t86 + t85 * t94;
	t62 = atan2(-t63, t74);
	t59 = sin(t62);
	t60 = cos(t62);
	t53 = -t59 * t63 + t60 * t74;
	t52 = 0.1e1 / t53 ^ 2;
	t73 = -t80 * t91 + t82 * t87;
	t66 = t73 * t84 - t80 * t93;
	t98 = t52 * t66;
	t67 = t73 * t86 + t80 * t94;
	t90 = t83 * t87;
	t72 = t80 * t90 + t82 * t85;
	t79 = qJ(4) + qJ(5);
	t77 = sin(t79);
	t78 = cos(t79);
	t58 = t67 * t78 + t72 * t77;
	t56 = 0.1e1 / t58 ^ 2;
	t57 = t67 * t77 - t72 * t78;
	t97 = t56 * t57;
	t69 = 0.1e1 / t74 ^ 2;
	t96 = t63 * t69;
	t95 = t72 * t86;
	t92 = t81 * t87;
	t89 = t57 ^ 2 * t56 + 0.1e1;
	t88 = -t59 * t74 - t60 * t63;
	t75 = t83 * t84 + t85 * t93;
	t70 = -t80 * t85 + t82 * t90;
	t68 = 0.1e1 / t74;
	t65 = t71 * t86 - t82 * t94;
	t61 = 0.1e1 / (t63 ^ 2 * t69 + 0.1e1);
	t55 = 0.1e1 / t58;
	t54 = 0.1e1 / t89;
	t51 = 0.1e1 / t53;
	t50 = 0.1e1 / (t66 ^ 2 * t52 + 0.1e1);
	t49 = (-t68 * t70 + t92 * t96) * t84 * t61;
	t48 = (-t65 * t68 + t75 * t96) * t61;
	t47 = t89 * t54;
	t1 = [0, t49, t48, 0, 0, 0; 0, (-t72 * t84 * t51 - ((-t59 * t70 + t60 * t92) * t84 + t88 * t49) * t98) * t50, (t67 * t51 - (t88 * t48 - t59 * t65 + t60 * t75) * t98) * t50, 0, 0, 0; 0, ((-t73 * t78 - t77 * t95) * t55 - (t73 * t77 - t78 * t95) * t97) * t54, (-t55 * t77 + t78 * t97) * t66 * t54, t47, t47, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:05:49
	% EndTime: 2019-10-09 23:05:49
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (427->31), mult. (1032->79), div. (70->9), fcn. (1461->13), ass. (0->51)
	t88 = sin(pkin(6));
	t93 = cos(qJ(3));
	t100 = t88 * t93;
	t87 = sin(pkin(11));
	t89 = cos(pkin(11));
	t94 = cos(qJ(2));
	t90 = cos(pkin(6));
	t92 = sin(qJ(2));
	t98 = t90 * t92;
	t78 = t87 * t94 + t89 * t98;
	t91 = sin(qJ(3));
	t70 = t89 * t100 + t78 * t91;
	t101 = t88 * t91;
	t81 = t92 * t101 - t90 * t93;
	t69 = atan2(-t70, t81);
	t66 = sin(t69);
	t67 = cos(t69);
	t60 = -t66 * t70 + t67 * t81;
	t59 = 0.1e1 / t60 ^ 2;
	t80 = -t87 * t98 + t89 * t94;
	t73 = -t87 * t100 + t80 * t91;
	t105 = t59 * t73;
	t74 = t87 * t101 + t80 * t93;
	t97 = t90 * t94;
	t79 = t87 * t97 + t89 * t92;
	t86 = qJ(4) + qJ(5);
	t84 = sin(t86);
	t85 = cos(t86);
	t65 = t74 * t85 + t79 * t84;
	t63 = 0.1e1 / t65 ^ 2;
	t64 = t74 * t84 - t79 * t85;
	t104 = t63 * t64;
	t76 = 0.1e1 / t81 ^ 2;
	t103 = t70 * t76;
	t102 = t79 * t93;
	t99 = t88 * t94;
	t96 = t64 ^ 2 * t63 + 0.1e1;
	t95 = -t66 * t81 - t67 * t70;
	t82 = t92 * t100 + t90 * t91;
	t77 = -t87 * t92 + t89 * t97;
	t75 = 0.1e1 / t81;
	t72 = -t89 * t101 + t78 * t93;
	t68 = 0.1e1 / (t70 ^ 2 * t76 + 0.1e1);
	t62 = 0.1e1 / t65;
	t61 = 0.1e1 / t96;
	t58 = 0.1e1 / t60;
	t57 = 0.1e1 / (t73 ^ 2 * t59 + 0.1e1);
	t56 = (t99 * t103 - t75 * t77) * t91 * t68;
	t55 = (t82 * t103 - t72 * t75) * t68;
	t54 = t96 * t61;
	t1 = [0, t56, t55, 0, 0, 0; 0, (-t79 * t91 * t58 - ((-t66 * t77 + t67 * t99) * t91 + t95 * t56) * t105) * t57, (t74 * t58 - (t95 * t55 - t66 * t72 + t67 * t82) * t105) * t57, 0, 0, 0; 0, ((-t84 * t102 - t80 * t85) * t62 - (-t85 * t102 + t80 * t84) * t104) * t61, (t85 * t104 - t62 * t84) * t73 * t61, t54, t54, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end