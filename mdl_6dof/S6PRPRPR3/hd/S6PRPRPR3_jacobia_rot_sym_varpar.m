% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR3
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
%   Wie in S6PRPRPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRPR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (30->0), mult. (78->0), div. (6->0), fcn. (108->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.17s
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
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (481->27), mult. (1334->73), div. (57->9), fcn. (1884->13), ass. (0->47)
	t78 = cos(pkin(6));
	t73 = sin(pkin(11));
	t76 = cos(pkin(11));
	t80 = sin(qJ(2));
	t82 = cos(qJ(2));
	t84 = t82 * t73 + t80 * t76;
	t69 = t84 * t78;
	t70 = t80 * t73 - t82 * t76;
	t74 = sin(pkin(10));
	t77 = cos(pkin(10));
	t58 = t77 * t69 - t74 * t70;
	t79 = sin(qJ(4));
	t75 = sin(pkin(6));
	t81 = cos(qJ(4));
	t86 = t75 * t81;
	t50 = t58 * t79 + t77 * t86;
	t68 = t84 * t75;
	t64 = t68 * t79 - t78 * t81;
	t49 = atan2(-t50, t64);
	t46 = sin(t49);
	t47 = cos(t49);
	t44 = -t46 * t50 + t47 * t64;
	t43 = 0.1e1 / t44 ^ 2;
	t61 = -t74 * t69 - t77 * t70;
	t53 = t61 * t79 - t74 * t86;
	t89 = t43 * t53;
	t63 = 0.1e1 / t64 ^ 2;
	t88 = t50 * t63;
	t87 = t75 * t79;
	t85 = -t46 * t64 - t47 * t50;
	t83 = t70 * t78;
	t67 = t70 * t75;
	t65 = t68 * t81 + t78 * t79;
	t62 = 0.1e1 / t64;
	t59 = t74 * t83 - t77 * t84;
	t57 = -t74 * t84 - t77 * t83;
	t56 = 0.1e1 / t59 ^ 2;
	t55 = 0.1e1 / t59;
	t54 = t61 * t81 + t74 * t87;
	t52 = t58 * t81 - t77 * t87;
	t48 = 0.1e1 / (t50 ^ 2 * t63 + 0.1e1);
	t45 = 0.1e1 / (t54 ^ 2 * t56 + 0.1e1);
	t42 = 0.1e1 / t44;
	t41 = 0.1e1 / (t53 ^ 2 * t43 + 0.1e1);
	t40 = (-t57 * t62 - t67 * t88) * t79 * t48;
	t39 = (-t52 * t62 + t65 * t88) * t48;
	t1 = [0, t40, 0, t39, 0, 0; 0, (t59 * t79 * t42 - ((-t46 * t57 - t47 * t67) * t79 + t85 * t40) * t89) * t41, 0, (t54 * t42 - (t85 * t39 - t46 * t52 + t47 * t65) * t89) * t41, 0, 0; 0, (-t61 * t54 * t56 - t59 * t81 * t55) * t45, 0, t53 * t55 * t45, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:33:40
	% EndTime: 2019-10-09 21:33:40
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (632->32), mult. (1727->83), div. (65->9), fcn. (2425->15), ass. (0->55)
	t87 = sin(pkin(6));
	t92 = sin(qJ(4));
	t103 = t87 * t92;
	t90 = cos(pkin(6));
	t85 = sin(pkin(11));
	t88 = cos(pkin(11));
	t93 = sin(qJ(2));
	t96 = cos(qJ(2));
	t98 = t96 * t85 + t93 * t88;
	t80 = t98 * t90;
	t81 = t93 * t85 - t96 * t88;
	t86 = sin(pkin(10));
	t89 = cos(pkin(10));
	t69 = t89 * t80 - t86 * t81;
	t95 = cos(qJ(4));
	t64 = -t89 * t103 + t69 * t95;
	t79 = t98 * t87;
	t76 = t79 * t95 + t90 * t92;
	t62 = atan2(-t64, t76);
	t59 = sin(t62);
	t60 = cos(t62);
	t53 = -t59 * t64 + t60 * t76;
	t52 = 0.1e1 / t53 ^ 2;
	t99 = -t86 * t80 - t89 * t81;
	t67 = t86 * t103 + t95 * t99;
	t108 = t52 * t67;
	t97 = t81 * t90;
	t71 = t86 * t97 - t89 * t98;
	t94 = cos(qJ(6));
	t104 = t71 * t94;
	t102 = t87 * t95;
	t66 = -t86 * t102 + t92 * t99;
	t91 = sin(qJ(6));
	t58 = t66 * t91 - t104;
	t56 = 0.1e1 / t58 ^ 2;
	t105 = t71 * t91;
	t57 = -t66 * t94 - t105;
	t107 = t56 * t57;
	t74 = 0.1e1 / t76 ^ 2;
	t106 = t64 * t74;
	t101 = t57 ^ 2 * t56 + 0.1e1;
	t100 = -t59 * t76 - t60 * t64;
	t78 = t81 * t87;
	t75 = -t79 * t92 + t90 * t95;
	t73 = 0.1e1 / t76;
	t68 = -t86 * t98 - t89 * t97;
	t63 = t89 * t102 + t69 * t92;
	t61 = 0.1e1 / (t64 ^ 2 * t74 + 0.1e1);
	t55 = 0.1e1 / t58;
	t54 = 0.1e1 / t101;
	t51 = 0.1e1 / t53;
	t50 = 0.1e1 / (t67 ^ 2 * t52 + 0.1e1);
	t49 = (-t78 * t106 - t68 * t73) * t95 * t61;
	t48 = (t75 * t106 + t63 * t73) * t61;
	t1 = [0, t49, 0, t48, 0, 0; 0, (t71 * t95 * t51 - ((-t59 * t68 - t60 * t78) * t95 + t100 * t49) * t108) * t50, 0, (-t66 * t51 - (t100 * t48 + t59 * t63 + t60 * t75) * t108) * t50, 0, 0; 0, ((-t92 * t104 + t91 * t99) * t55 - (t92 * t105 + t94 * t99) * t107) * t54, 0, (-t91 * t107 - t55 * t94) * t67 * t54, 0, t101 * t54;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end