% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR1
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
%   Wie in S6PRPRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->0), mult. (78->0), div. (6->0), fcn. (108->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (222->18), mult. (626->46), div. (35->9), fcn. (880->13), ass. (0->34)
	t56 = sin(pkin(11));
	t57 = sin(pkin(6));
	t68 = t56 * t57;
	t60 = cos(pkin(6));
	t55 = sin(pkin(12));
	t58 = cos(pkin(12));
	t62 = sin(qJ(2));
	t64 = cos(qJ(2));
	t66 = t64 * t55 + t62 * t58;
	t52 = t66 * t60;
	t53 = t62 * t55 - t64 * t58;
	t59 = cos(pkin(11));
	t47 = -t56 * t52 - t59 * t53;
	t61 = sin(qJ(4));
	t63 = cos(qJ(4));
	t42 = t47 * t63 + t61 * t68;
	t40 = 0.1e1 / t42 ^ 2;
	t41 = t47 * t61 - t63 * t68;
	t67 = t41 ^ 2 * t40 + 0.1e1;
	t65 = t53 * t60;
	t51 = t66 * t57;
	t50 = t53 * t57;
	t49 = 0.1e1 / t50 ^ 2;
	t45 = t56 * t65 - t59 * t66;
	t44 = -t56 * t66 - t59 * t65;
	t43 = -t59 * t52 + t56 * t53;
	t39 = atan2(t44, t50);
	t37 = cos(t39);
	t36 = sin(t39);
	t35 = 0.1e1 / t67;
	t34 = t36 * t44 + t37 * t50;
	t33 = 0.1e1 / t34 ^ 2;
	t31 = (t43 / t50 - t51 * t44 * t49) / (t44 ^ 2 * t49 + 0.1e1);
	t1 = [0, t31, 0, 0, 0, 0; 0, (t47 / t34 + (t36 * t43 + t37 * t51 + (-t36 * t50 + t37 * t44) * t31) * t45 * t33) / (t45 ^ 2 * t33 + 0.1e1), 0, 0, 0, 0; 0, (t61 / t42 - t63 * t41 * t40) * t45 * t35, 0, t67 * t35, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (297->19), mult. (709->46), div. (40->9), fcn. (992->13), ass. (0->36)
	t74 = sin(pkin(11));
	t75 = sin(pkin(6));
	t84 = t74 * t75;
	t78 = cos(pkin(6));
	t73 = sin(pkin(12));
	t76 = cos(pkin(12));
	t79 = sin(qJ(2));
	t80 = cos(qJ(2));
	t82 = t80 * t73 + t79 * t76;
	t67 = t82 * t78;
	t68 = t79 * t73 - t80 * t76;
	t77 = cos(pkin(11));
	t62 = -t74 * t67 - t77 * t68;
	t72 = qJ(4) + qJ(5);
	t70 = sin(t72);
	t71 = cos(t72);
	t57 = t62 * t71 + t70 * t84;
	t55 = 0.1e1 / t57 ^ 2;
	t56 = t62 * t70 - t71 * t84;
	t83 = t56 ^ 2 * t55 + 0.1e1;
	t81 = t68 * t78;
	t66 = t82 * t75;
	t65 = t68 * t75;
	t64 = 0.1e1 / t65 ^ 2;
	t60 = t74 * t81 - t77 * t82;
	t59 = -t74 * t82 - t77 * t81;
	t58 = -t77 * t67 + t74 * t68;
	t54 = atan2(t59, t65);
	t52 = cos(t54);
	t51 = sin(t54);
	t50 = 0.1e1 / t83;
	t49 = t51 * t59 + t52 * t65;
	t48 = 0.1e1 / t49 ^ 2;
	t46 = (t58 / t65 - t66 * t59 * t64) / (t59 ^ 2 * t64 + 0.1e1);
	t45 = t83 * t50;
	t1 = [0, t46, 0, 0, 0, 0; 0, (t62 / t49 + (t51 * t58 + t52 * t66 + (-t51 * t65 + t52 * t59) * t46) * t60 * t48) / (t60 ^ 2 * t48 + 0.1e1), 0, 0, 0, 0; 0, (t70 / t57 - t71 * t56 * t55) * t60 * t50, 0, t45, t45, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (1379->33), mult. (2514->83), div. (95->9), fcn. (3534->15), ass. (0->58)
	t109 = qJ(4) + qJ(5);
	t107 = sin(t109);
	t108 = cos(t109);
	t112 = sin(pkin(6));
	t114 = cos(pkin(11));
	t125 = t112 * t114;
	t115 = cos(pkin(6));
	t110 = sin(pkin(12));
	t113 = cos(pkin(12));
	t117 = sin(qJ(2));
	t119 = cos(qJ(2));
	t121 = t119 * t110 + t117 * t113;
	t102 = t121 * t115;
	t103 = t117 * t110 - t119 * t113;
	t111 = sin(pkin(11));
	t91 = t114 * t102 - t111 * t103;
	t85 = t91 * t107 + t108 * t125;
	t101 = t121 * t112;
	t97 = t101 * t107 - t115 * t108;
	t84 = atan2(-t85, t97);
	t81 = sin(t84);
	t82 = cos(t84);
	t75 = -t81 * t85 + t82 * t97;
	t74 = 0.1e1 / t75 ^ 2;
	t122 = -t111 * t102 - t114 * t103;
	t126 = t111 * t112;
	t88 = t107 * t122 - t108 * t126;
	t131 = t74 * t88;
	t118 = cos(qJ(6));
	t116 = sin(qJ(6));
	t120 = t103 * t115;
	t93 = t111 * t120 - t114 * t121;
	t128 = t93 * t116;
	t89 = t107 * t126 + t108 * t122;
	t80 = t89 * t118 - t128;
	t78 = 0.1e1 / t80 ^ 2;
	t127 = t93 * t118;
	t79 = t89 * t116 + t127;
	t130 = t78 * t79;
	t96 = 0.1e1 / t97 ^ 2;
	t129 = t85 * t96;
	t124 = t79 ^ 2 * t78 + 0.1e1;
	t123 = -t81 * t97 - t82 * t85;
	t100 = t103 * t112;
	t98 = t101 * t108 + t115 * t107;
	t95 = 0.1e1 / t97;
	t90 = -t111 * t121 - t114 * t120;
	t87 = -t107 * t125 + t91 * t108;
	t83 = 0.1e1 / (t85 ^ 2 * t96 + 0.1e1);
	t77 = 0.1e1 / t80;
	t76 = 0.1e1 / t124;
	t73 = 0.1e1 / t75;
	t72 = 0.1e1 / (t88 ^ 2 * t74 + 0.1e1);
	t71 = (-t100 * t129 - t90 * t95) * t83 * t107;
	t70 = (t98 * t129 - t87 * t95) * t83;
	t69 = (-t116 * t77 + t118 * t130) * t88 * t76;
	t68 = (t89 * t73 - (t123 * t70 - t81 * t87 + t82 * t98) * t131) * t72;
	t1 = [0, t71, 0, t70, t70, 0; 0, (t93 * t107 * t73 - (t123 * t71 + (-t100 * t82 - t81 * t90) * t107) * t131) * t72, 0, t68, t68, 0; 0, ((t108 * t128 - t118 * t122) * t77 - (t108 * t127 + t116 * t122) * t130) * t76, 0, t69, t69, t124 * t76;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end