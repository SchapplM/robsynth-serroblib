% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR1
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
%   Wie in S6PRPRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->0), mult. (78->0), div. (6->0), fcn. (108->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (222->18), mult. (626->46), div. (35->9), fcn. (880->13), ass. (0->34)
	t56 = sin(pkin(10));
	t57 = sin(pkin(6));
	t68 = t56 * t57;
	t60 = cos(pkin(6));
	t55 = sin(pkin(11));
	t58 = cos(pkin(11));
	t62 = sin(qJ(2));
	t64 = cos(qJ(2));
	t66 = t64 * t55 + t62 * t58;
	t52 = t66 * t60;
	t53 = t62 * t55 - t64 * t58;
	t59 = cos(pkin(10));
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
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (252->19), mult. (626->46), div. (35->9), fcn. (880->13), ass. (0->35)
	t67 = sin(pkin(10));
	t68 = sin(pkin(6));
	t77 = t67 * t68;
	t71 = cos(pkin(6));
	t66 = sin(pkin(11));
	t69 = cos(pkin(11));
	t72 = sin(qJ(2));
	t73 = cos(qJ(2));
	t75 = t73 * t66 + t72 * t69;
	t60 = t75 * t71;
	t61 = t72 * t66 - t73 * t69;
	t70 = cos(pkin(10));
	t55 = -t67 * t60 - t70 * t61;
	t65 = qJ(4) + pkin(12);
	t63 = sin(t65);
	t64 = cos(t65);
	t50 = t55 * t64 + t63 * t77;
	t48 = 0.1e1 / t50 ^ 2;
	t49 = t55 * t63 - t64 * t77;
	t76 = t49 ^ 2 * t48 + 0.1e1;
	t74 = t61 * t71;
	t59 = t75 * t68;
	t58 = t61 * t68;
	t57 = 0.1e1 / t58 ^ 2;
	t53 = t67 * t74 - t70 * t75;
	t52 = -t67 * t75 - t70 * t74;
	t51 = -t70 * t60 + t67 * t61;
	t47 = atan2(t52, t58);
	t45 = cos(t47);
	t44 = sin(t47);
	t43 = 0.1e1 / t76;
	t42 = t44 * t52 + t45 * t58;
	t41 = 0.1e1 / t42 ^ 2;
	t39 = (t51 / t58 - t59 * t52 * t57) / (t52 ^ 2 * t57 + 0.1e1);
	t1 = [0, t39, 0, 0, 0, 0; 0, (t55 / t42 + (t44 * t51 + t45 * t59 + (-t44 * t58 + t45 * t52) * t39) * t53 * t41) / (t53 ^ 2 * t41 + 0.1e1), 0, 0, 0, 0; 0, (t63 / t50 - t64 * t49 * t48) * t53 * t43, 0, t76 * t43, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:29:59
	% EndTime: 2019-10-09 21:29:59
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (939->33), mult. (1727->83), div. (65->9), fcn. (2425->15), ass. (0->56)
	t93 = sin(pkin(6));
	t95 = cos(pkin(10));
	t107 = t93 * t95;
	t100 = cos(qJ(2));
	t91 = sin(pkin(11));
	t94 = cos(pkin(11));
	t98 = sin(qJ(2));
	t103 = t100 * t94 - t98 * t91;
	t102 = t100 * t91 + t98 * t94;
	t96 = cos(pkin(6));
	t83 = t102 * t96;
	t92 = sin(pkin(10));
	t72 = t103 * t92 + t95 * t83;
	t90 = qJ(4) + pkin(12);
	t88 = sin(t90);
	t89 = cos(t90);
	t66 = t89 * t107 + t72 * t88;
	t82 = t102 * t93;
	t78 = t82 * t88 - t96 * t89;
	t65 = atan2(-t66, t78);
	t62 = sin(t65);
	t63 = cos(t65);
	t56 = -t62 * t66 + t63 * t78;
	t55 = 0.1e1 / t56 ^ 2;
	t104 = t103 * t95 - t92 * t83;
	t108 = t92 * t93;
	t69 = t104 * t88 - t89 * t108;
	t113 = t55 * t69;
	t101 = t103 * t96;
	t74 = -t92 * t101 - t102 * t95;
	t97 = sin(qJ(6));
	t110 = t74 * t97;
	t70 = t104 * t89 + t88 * t108;
	t99 = cos(qJ(6));
	t61 = t70 * t99 - t110;
	t59 = 0.1e1 / t61 ^ 2;
	t109 = t74 * t99;
	t60 = t70 * t97 + t109;
	t112 = t59 * t60;
	t77 = 0.1e1 / t78 ^ 2;
	t111 = t66 * t77;
	t106 = t60 ^ 2 * t59 + 0.1e1;
	t105 = -t62 * t78 - t63 * t66;
	t81 = t103 * t93;
	t79 = t82 * t89 + t96 * t88;
	t76 = 0.1e1 / t78;
	t71 = t95 * t101 - t102 * t92;
	t68 = -t88 * t107 + t72 * t89;
	t64 = 0.1e1 / (t66 ^ 2 * t77 + 0.1e1);
	t58 = 0.1e1 / t61;
	t57 = 0.1e1 / t106;
	t54 = 0.1e1 / t56;
	t53 = 0.1e1 / (t69 ^ 2 * t55 + 0.1e1);
	t52 = (t81 * t111 - t71 * t76) * t88 * t64;
	t51 = (t79 * t111 - t68 * t76) * t64;
	t1 = [0, t52, 0, t51, 0, 0; 0, (t74 * t88 * t54 - ((-t62 * t71 + t63 * t81) * t88 + t105 * t52) * t113) * t53, 0, (t70 * t54 - (t105 * t51 - t62 * t68 + t63 * t79) * t113) * t53, 0, 0; 0, ((-t104 * t99 + t89 * t110) * t58 - (t104 * t97 + t89 * t109) * t112) * t57, 0, (t99 * t112 - t58 * t97) * t69 * t57, 0, t106 * t57;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end