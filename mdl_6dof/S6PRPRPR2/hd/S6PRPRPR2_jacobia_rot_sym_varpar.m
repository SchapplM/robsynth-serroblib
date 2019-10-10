% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR2
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
%   Wie in S6PRPRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->0), mult. (78->0), div. (6->0), fcn. (108->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
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
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (573->32), mult. (1574->83), div. (60->9), fcn. (2215->15), ass. (0->53)
	t88 = cos(pkin(6));
	t82 = sin(pkin(11));
	t86 = cos(pkin(11));
	t90 = sin(qJ(2));
	t92 = cos(qJ(2));
	t94 = t82 * t92 + t86 * t90;
	t76 = t94 * t88;
	t77 = t82 * t90 - t92 * t86;
	t83 = sin(pkin(10));
	t87 = cos(pkin(10));
	t65 = t76 * t87 - t77 * t83;
	t89 = sin(qJ(4));
	t84 = sin(pkin(6));
	t91 = cos(qJ(4));
	t97 = t84 * t91;
	t59 = t65 * t89 + t87 * t97;
	t75 = t94 * t84;
	t71 = t75 * t89 - t88 * t91;
	t58 = atan2(-t59, t71);
	t55 = sin(t58);
	t56 = cos(t58);
	t49 = -t55 * t59 + t56 * t71;
	t48 = 0.1e1 / t49 ^ 2;
	t95 = -t76 * t83 - t77 * t87;
	t62 = -t83 * t97 + t89 * t95;
	t102 = t48 * t62;
	t98 = t84 * t89;
	t63 = t83 * t98 + t91 * t95;
	t93 = t77 * t88;
	t67 = t83 * t93 - t87 * t94;
	t81 = sin(pkin(12));
	t85 = cos(pkin(12));
	t54 = t63 * t85 - t67 * t81;
	t52 = 0.1e1 / t54 ^ 2;
	t53 = t63 * t81 + t67 * t85;
	t101 = t52 * t53;
	t70 = 0.1e1 / t71 ^ 2;
	t100 = t59 * t70;
	t99 = t67 * t91;
	t96 = -t55 * t71 - t56 * t59;
	t74 = t77 * t84;
	t72 = t75 * t91 + t88 * t89;
	t69 = 0.1e1 / t71;
	t64 = -t83 * t94 - t87 * t93;
	t61 = t65 * t91 - t87 * t98;
	t57 = 0.1e1 / (t59 ^ 2 * t70 + 0.1e1);
	t51 = 0.1e1 / t54;
	t50 = 0.1e1 / (t52 * t53 ^ 2 + 0.1e1);
	t47 = 0.1e1 / t49;
	t46 = 0.1e1 / (t48 * t62 ^ 2 + 0.1e1);
	t45 = (-t100 * t74 - t64 * t69) * t89 * t57;
	t44 = (t100 * t72 - t61 * t69) * t57;
	t1 = [0, t45, 0, t44, 0, 0; 0, (t67 * t89 * t47 - ((-t55 * t64 - t56 * t74) * t89 + t96 * t45) * t102) * t46, 0, (t63 * t47 - (t44 * t96 - t55 * t61 + t56 * t72) * t102) * t46, 0, 0; 0, ((t81 * t99 - t85 * t95) * t51 - (t81 * t95 + t85 * t99) * t101) * t50, 0, (t101 * t85 - t51 * t81) * t62 * t50, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (680->33), mult. (1727->84), div. (65->9), fcn. (2425->15), ass. (0->55)
	t100 = cos(qJ(4));
	t94 = sin(pkin(6));
	t108 = t100 * t94;
	t101 = cos(qJ(2));
	t92 = sin(pkin(11));
	t95 = cos(pkin(11));
	t99 = sin(qJ(2));
	t104 = t101 * t95 - t99 * t92;
	t103 = t101 * t92 + t99 * t95;
	t97 = cos(pkin(6));
	t84 = t103 * t97;
	t93 = sin(pkin(10));
	t96 = cos(pkin(10));
	t73 = t104 * t93 + t96 * t84;
	t98 = sin(qJ(4));
	t67 = t96 * t108 + t73 * t98;
	t83 = t103 * t94;
	t79 = -t97 * t100 + t83 * t98;
	t66 = atan2(-t67, t79);
	t63 = sin(t66);
	t64 = cos(t66);
	t57 = -t63 * t67 + t64 * t79;
	t56 = 0.1e1 / t57 ^ 2;
	t105 = t104 * t96 - t93 * t84;
	t70 = t105 * t98 - t93 * t108;
	t113 = t56 * t70;
	t110 = t94 * t98;
	t71 = t100 * t105 + t93 * t110;
	t102 = t104 * t97;
	t75 = -t93 * t102 - t103 * t96;
	t91 = pkin(12) + qJ(6);
	t89 = sin(t91);
	t90 = cos(t91);
	t62 = t71 * t90 - t75 * t89;
	t60 = 0.1e1 / t62 ^ 2;
	t61 = t71 * t89 + t75 * t90;
	t112 = t60 * t61;
	t78 = 0.1e1 / t79 ^ 2;
	t111 = t67 * t78;
	t109 = t100 * t75;
	t107 = t61 ^ 2 * t60 + 0.1e1;
	t106 = -t63 * t79 - t64 * t67;
	t82 = t104 * t94;
	t80 = t83 * t100 + t97 * t98;
	t77 = 0.1e1 / t79;
	t72 = t96 * t102 - t103 * t93;
	t69 = t73 * t100 - t96 * t110;
	t65 = 0.1e1 / (t67 ^ 2 * t78 + 0.1e1);
	t59 = 0.1e1 / t62;
	t58 = 0.1e1 / t107;
	t55 = 0.1e1 / t57;
	t54 = 0.1e1 / (t70 ^ 2 * t56 + 0.1e1);
	t53 = (t82 * t111 - t72 * t77) * t98 * t65;
	t52 = (t80 * t111 - t69 * t77) * t65;
	t1 = [0, t53, 0, t52, 0, 0; 0, (t75 * t98 * t55 - ((-t63 * t72 + t64 * t82) * t98 + t106 * t53) * t113) * t54, 0, (t71 * t55 - (t106 * t52 - t63 * t69 + t64 * t80) * t113) * t54, 0, 0; 0, ((-t105 * t90 + t89 * t109) * t59 - (t105 * t89 + t90 * t109) * t112) * t58, 0, (t90 * t112 - t59 * t89) * t70 * t58, 0, t107 * t58;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end