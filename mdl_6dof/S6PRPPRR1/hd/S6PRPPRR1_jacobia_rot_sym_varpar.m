% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR1
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
%   Wie in S6PRPPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:24
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPPRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobia_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->0), mult. (78->0), div. (6->0), fcn. (108->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (191->18), mult. (543->44), div. (30->9), fcn. (768->13), ass. (0->32)
	t54 = sin(pkin(10));
	t55 = sin(pkin(6));
	t64 = t54 * t55;
	t59 = cos(pkin(6));
	t53 = sin(pkin(11));
	t57 = cos(pkin(11));
	t60 = sin(qJ(2));
	t61 = cos(qJ(2));
	t63 = t61 * t53 + t60 * t57;
	t49 = t63 * t59;
	t50 = t60 * t53 - t61 * t57;
	t58 = cos(pkin(10));
	t44 = -t54 * t49 - t58 * t50;
	t62 = t50 * t59;
	t56 = cos(pkin(12));
	t52 = sin(pkin(12));
	t48 = t63 * t55;
	t47 = t50 * t55;
	t46 = 0.1e1 / t47 ^ 2;
	t42 = t54 * t62 - t58 * t63;
	t41 = -t54 * t63 - t58 * t62;
	t40 = -t58 * t49 + t54 * t50;
	t39 = t44 * t56 + t52 * t64;
	t38 = t44 * t52 - t56 * t64;
	t37 = 0.1e1 / t39 ^ 2;
	t36 = atan2(t41, t47);
	t34 = cos(t36);
	t33 = sin(t36);
	t31 = t33 * t41 + t34 * t47;
	t30 = 0.1e1 / t31 ^ 2;
	t28 = (t40 / t47 - t48 * t41 * t46) / (t41 ^ 2 * t46 + 0.1e1);
	t1 = [0, t28, 0, 0, 0, 0; 0, (t44 / t31 + (t33 * t40 + t34 * t48 + (-t33 * t47 + t34 * t41) * t28) * t42 * t30) / (t42 ^ 2 * t30 + 0.1e1), 0, 0, 0, 0; 0, (t52 / t39 - t56 * t38 * t37) * t42 / (t38 ^ 2 * t37 + 0.1e1), 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (252->19), mult. (626->46), div. (35->9), fcn. (880->13), ass. (0->35)
	t65 = sin(pkin(10));
	t66 = sin(pkin(6));
	t75 = t65 * t66;
	t69 = cos(pkin(6));
	t64 = sin(pkin(11));
	t67 = cos(pkin(11));
	t70 = sin(qJ(2));
	t71 = cos(qJ(2));
	t73 = t71 * t64 + t70 * t67;
	t58 = t73 * t69;
	t59 = t70 * t64 - t71 * t67;
	t68 = cos(pkin(10));
	t53 = -t65 * t58 - t68 * t59;
	t63 = pkin(12) + qJ(5);
	t61 = sin(t63);
	t62 = cos(t63);
	t48 = t53 * t62 + t61 * t75;
	t46 = 0.1e1 / t48 ^ 2;
	t47 = t53 * t61 - t62 * t75;
	t74 = t47 ^ 2 * t46 + 0.1e1;
	t72 = t59 * t69;
	t57 = t73 * t66;
	t56 = t59 * t66;
	t55 = 0.1e1 / t56 ^ 2;
	t51 = t65 * t72 - t68 * t73;
	t50 = -t65 * t73 - t68 * t72;
	t49 = -t68 * t58 + t65 * t59;
	t45 = atan2(t50, t56);
	t43 = cos(t45);
	t42 = sin(t45);
	t41 = 0.1e1 / t74;
	t40 = t42 * t50 + t43 * t56;
	t39 = 0.1e1 / t40 ^ 2;
	t37 = (t49 / t56 - t57 * t50 * t55) / (t50 ^ 2 * t55 + 0.1e1);
	t1 = [0, t37, 0, 0, 0, 0; 0, (t53 / t40 + (t42 * t49 + t43 * t57 + (-t42 * t56 + t43 * t50) * t37) * t51 * t39) / (t51 ^ 2 * t39 + 0.1e1), 0, 0, 0, 0; 0, (t61 / t48 - t62 * t47 * t46) * t51 * t41, 0, 0, t74 * t41, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (939->33), mult. (1727->83), div. (65->9), fcn. (2425->15), ass. (0->56)
	t91 = sin(pkin(6));
	t93 = cos(pkin(10));
	t104 = t91 * t93;
	t89 = sin(pkin(11));
	t92 = cos(pkin(11));
	t96 = sin(qJ(2));
	t98 = cos(qJ(2));
	t100 = t98 * t89 + t96 * t92;
	t94 = cos(pkin(6));
	t81 = t100 * t94;
	t82 = t96 * t89 - t98 * t92;
	t90 = sin(pkin(10));
	t70 = t93 * t81 - t90 * t82;
	t88 = pkin(12) + qJ(5);
	t86 = sin(t88);
	t87 = cos(t88);
	t64 = t87 * t104 + t70 * t86;
	t80 = t100 * t91;
	t76 = t80 * t86 - t94 * t87;
	t63 = atan2(-t64, t76);
	t60 = sin(t63);
	t61 = cos(t63);
	t54 = -t60 * t64 + t61 * t76;
	t53 = 0.1e1 / t54 ^ 2;
	t101 = -t90 * t81 - t93 * t82;
	t105 = t90 * t91;
	t67 = t101 * t86 - t87 * t105;
	t110 = t53 * t67;
	t99 = t82 * t94;
	t72 = -t100 * t93 + t90 * t99;
	t95 = sin(qJ(6));
	t107 = t72 * t95;
	t68 = t101 * t87 + t86 * t105;
	t97 = cos(qJ(6));
	t59 = t68 * t97 - t107;
	t57 = 0.1e1 / t59 ^ 2;
	t106 = t72 * t97;
	t58 = t68 * t95 + t106;
	t109 = t57 * t58;
	t75 = 0.1e1 / t76 ^ 2;
	t108 = t64 * t75;
	t103 = t58 ^ 2 * t57 + 0.1e1;
	t102 = -t60 * t76 - t61 * t64;
	t79 = t82 * t91;
	t77 = t80 * t87 + t94 * t86;
	t74 = 0.1e1 / t76;
	t69 = -t100 * t90 - t93 * t99;
	t66 = -t86 * t104 + t70 * t87;
	t62 = 0.1e1 / (t64 ^ 2 * t75 + 0.1e1);
	t56 = 0.1e1 / t59;
	t55 = 0.1e1 / t103;
	t52 = 0.1e1 / t54;
	t51 = 0.1e1 / (t67 ^ 2 * t53 + 0.1e1);
	t50 = (-t79 * t108 - t69 * t74) * t86 * t62;
	t49 = (t77 * t108 - t66 * t74) * t62;
	t1 = [0, t50, 0, 0, t49, 0; 0, (t72 * t86 * t52 - ((-t60 * t69 - t61 * t79) * t86 + t102 * t50) * t110) * t51, 0, 0, (t68 * t52 - (t102 * t49 - t60 * t66 + t61 * t77) * t110) * t51, 0; 0, ((-t101 * t97 + t87 * t107) * t56 - (t101 * t95 + t87 * t106) * t109) * t55, 0, 0, (t97 * t109 - t56 * t95) * t67 * t55, t103 * t55;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end