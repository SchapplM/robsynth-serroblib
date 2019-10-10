% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR4
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
%   Wie in S6RRRRPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:36
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRPR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:50
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:50
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:50
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (358->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t68 = qJ(2) + qJ(3);
	t67 = cos(t68);
	t66 = sin(t68);
	t70 = sin(qJ(1));
	t78 = t70 * t66;
	t61 = atan2(-t78, -t67);
	t59 = sin(t61);
	t60 = cos(t61);
	t52 = -t59 * t78 - t60 * t67;
	t51 = 0.1e1 / t52 ^ 2;
	t72 = cos(qJ(1));
	t84 = t51 * t72 ^ 2;
	t71 = cos(qJ(4));
	t74 = t72 * t71;
	t69 = sin(qJ(4));
	t77 = t70 * t69;
	t58 = t67 * t74 + t77;
	t56 = 0.1e1 / t58 ^ 2;
	t75 = t72 * t69;
	t76 = t70 * t71;
	t57 = t67 * t75 - t76;
	t83 = t56 * t57;
	t82 = t59 * t67;
	t63 = t66 ^ 2;
	t81 = t63 / t67 ^ 2;
	t80 = t66 * t72;
	t62 = 0.1e1 / (t70 ^ 2 * t81 + 0.1e1);
	t79 = t70 * t62;
	t73 = t57 ^ 2 * t56 + 0.1e1;
	t64 = 0.1e1 / t67;
	t55 = 0.1e1 / t58;
	t54 = 0.1e1 / t73;
	t53 = (0.1e1 + t81) * t79;
	t50 = 0.1e1 / t52;
	t49 = 0.1e1 / (t63 * t84 + 0.1e1);
	t48 = (-t55 * t69 + t71 * t83) * t54 * t80;
	t47 = (t67 * t50 - (-t70 * t82 + t60 * t66 + (-t60 * t78 + t82) * t53) * t66 * t51) * t72 * t49;
	t1 = [t64 * t62 * t80, t53, t53, 0, 0, 0; (-t50 * t78 - (-t60 * t63 * t64 * t79 + (t62 - 0.1e1) * t66 * t59) * t66 * t84) * t49, t47, t47, 0, 0, 0; ((-t67 * t77 - t74) * t55 - (-t67 * t76 + t75) * t83) * t54, t48, t48, t73 * t54, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (422->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
	t72 = qJ(2) + qJ(3);
	t70 = cos(t72);
	t69 = sin(t72);
	t73 = sin(qJ(1));
	t78 = t73 * t69;
	t62 = atan2(-t78, -t70);
	t60 = sin(t62);
	t61 = cos(t62);
	t53 = -t60 * t78 - t61 * t70;
	t52 = 0.1e1 / t53 ^ 2;
	t74 = cos(qJ(1));
	t86 = t52 * t74 ^ 2;
	t71 = qJ(4) + pkin(11);
	t68 = cos(t71);
	t76 = t74 * t68;
	t67 = sin(t71);
	t80 = t73 * t67;
	t59 = t70 * t76 + t80;
	t57 = 0.1e1 / t59 ^ 2;
	t77 = t74 * t67;
	t79 = t73 * t68;
	t58 = t70 * t77 - t79;
	t85 = t57 * t58;
	t84 = t60 * t70;
	t64 = t69 ^ 2;
	t83 = t64 / t70 ^ 2;
	t82 = t69 * t74;
	t63 = 0.1e1 / (t73 ^ 2 * t83 + 0.1e1);
	t81 = t73 * t63;
	t75 = t58 ^ 2 * t57 + 0.1e1;
	t65 = 0.1e1 / t70;
	t56 = 0.1e1 / t59;
	t55 = (0.1e1 + t83) * t81;
	t54 = 0.1e1 / t75;
	t51 = 0.1e1 / t53;
	t50 = 0.1e1 / (t64 * t86 + 0.1e1);
	t49 = (-t56 * t67 + t68 * t85) * t54 * t82;
	t48 = (t70 * t51 - (-t73 * t84 + t61 * t69 + (-t61 * t78 + t84) * t55) * t69 * t52) * t74 * t50;
	t1 = [t65 * t63 * t82, t55, t55, 0, 0, 0; (-t51 * t78 - (-t61 * t64 * t65 * t81 + (t63 - 0.1e1) * t69 * t60) * t69 * t86) * t50, t48, t48, 0, 0, 0; ((-t70 * t80 - t76) * t56 - (-t70 * t79 + t77) * t85) * t54, t49, t49, t75 * t54, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:36:51
	% EndTime: 2019-10-10 12:36:51
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (531->22), mult. (359->54), div. (84->9), fcn. (524->9), ass. (0->40)
	t82 = qJ(2) + qJ(3);
	t81 = cos(t82);
	t80 = sin(t82);
	t83 = sin(qJ(1));
	t88 = t83 * t80;
	t72 = atan2(-t88, -t81);
	t70 = sin(t72);
	t71 = cos(t72);
	t64 = -t70 * t88 - t71 * t81;
	t63 = 0.1e1 / t64 ^ 2;
	t84 = cos(qJ(1));
	t96 = t63 * t84 ^ 2;
	t79 = qJ(4) + pkin(11) + qJ(6);
	t75 = cos(t79);
	t86 = t84 * t75;
	t74 = sin(t79);
	t90 = t83 * t74;
	t69 = t81 * t86 + t90;
	t67 = 0.1e1 / t69 ^ 2;
	t87 = t84 * t74;
	t89 = t83 * t75;
	t68 = t81 * t87 - t89;
	t95 = t67 * t68;
	t94 = t70 * t81;
	t76 = t80 ^ 2;
	t93 = t76 / t81 ^ 2;
	t92 = t80 * t84;
	t73 = 0.1e1 / (t83 ^ 2 * t93 + 0.1e1);
	t91 = t83 * t73;
	t85 = t68 ^ 2 * t67 + 0.1e1;
	t77 = 0.1e1 / t81;
	t66 = 0.1e1 / t69;
	t65 = (0.1e1 + t93) * t91;
	t62 = 0.1e1 / t64;
	t61 = 0.1e1 / t85;
	t60 = 0.1e1 / (t76 * t96 + 0.1e1);
	t59 = t85 * t61;
	t58 = (-t66 * t74 + t75 * t95) * t61 * t92;
	t57 = (t81 * t62 - (-t83 * t94 + t71 * t80 + (-t71 * t88 + t94) * t65) * t80 * t63) * t84 * t60;
	t1 = [t77 * t73 * t92, t65, t65, 0, 0, 0; (-t62 * t88 - (-t71 * t76 * t77 * t91 + (t73 - 0.1e1) * t80 * t70) * t80 * t96) * t60, t57, t57, 0, 0, 0; ((-t81 * t90 - t86) * t66 - (-t81 * t89 + t87) * t95) * t61, t58, t58, t59, 0, t59;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end