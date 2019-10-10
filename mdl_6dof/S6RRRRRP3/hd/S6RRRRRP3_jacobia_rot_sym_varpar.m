% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRP3
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
%   Wie in S6RRRRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:00
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:00:21
	% EndTime: 2019-10-10 13:00:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:00:21
	% EndTime: 2019-10-10 13:00:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:00:21
	% EndTime: 2019-10-10 13:00:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:00:21
	% EndTime: 2019-10-10 13:00:21
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:00:21
	% EndTime: 2019-10-10 13:00:21
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
	% StartTime: 2019-10-10 13:00:21
	% EndTime: 2019-10-10 13:00:21
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (453->22), mult. (359->54), div. (84->9), fcn. (524->9), ass. (0->40)
	t80 = qJ(2) + qJ(3);
	t78 = cos(t80);
	t76 = sin(t80);
	t81 = sin(qJ(1));
	t87 = t81 * t76;
	t70 = atan2(-t87, -t78);
	t68 = sin(t70);
	t69 = cos(t70);
	t61 = -t68 * t87 - t69 * t78;
	t60 = 0.1e1 / t61 ^ 2;
	t82 = cos(qJ(1));
	t94 = t60 * t82 ^ 2;
	t79 = qJ(4) + qJ(5);
	t77 = cos(t79);
	t84 = t82 * t77;
	t75 = sin(t79);
	t88 = t81 * t75;
	t67 = t78 * t84 + t88;
	t65 = 0.1e1 / t67 ^ 2;
	t85 = t82 * t75;
	t86 = t81 * t77;
	t66 = t78 * t85 - t86;
	t93 = t65 * t66;
	t92 = t68 * t78;
	t72 = t76 ^ 2;
	t91 = t72 / t78 ^ 2;
	t90 = t76 * t82;
	t71 = 0.1e1 / (t81 ^ 2 * t91 + 0.1e1);
	t89 = t81 * t71;
	t83 = t66 ^ 2 * t65 + 0.1e1;
	t73 = 0.1e1 / t78;
	t64 = 0.1e1 / t67;
	t63 = (0.1e1 + t91) * t89;
	t62 = 0.1e1 / t83;
	t59 = 0.1e1 / t61;
	t58 = 0.1e1 / (t72 * t94 + 0.1e1);
	t57 = t83 * t62;
	t56 = (-t64 * t75 + t77 * t93) * t62 * t90;
	t55 = (t78 * t59 - (-t81 * t92 + t69 * t76 + (-t69 * t87 + t92) * t63) * t76 * t60) * t82 * t58;
	t1 = [t73 * t71 * t90, t63, t63, 0, 0, 0; (-t59 * t87 - (-t69 * t72 * t73 * t89 + (t71 - 0.1e1) * t76 * t68) * t76 * t94) * t58, t55, t55, 0, 0, 0; ((-t78 * t88 - t84) * t64 - (-t78 * t86 + t85) * t93) * t62, t56, t56, t57, t57, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:00:21
	% EndTime: 2019-10-10 13:00:21
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (453->22), mult. (359->54), div. (84->9), fcn. (524->9), ass. (0->40)
	t81 = qJ(2) + qJ(3);
	t79 = cos(t81);
	t77 = sin(t81);
	t82 = sin(qJ(1));
	t88 = t82 * t77;
	t71 = atan2(-t88, -t79);
	t69 = sin(t71);
	t70 = cos(t71);
	t62 = -t69 * t88 - t70 * t79;
	t61 = 0.1e1 / t62 ^ 2;
	t83 = cos(qJ(1));
	t95 = t61 * t83 ^ 2;
	t80 = qJ(4) + qJ(5);
	t78 = cos(t80);
	t85 = t83 * t78;
	t76 = sin(t80);
	t89 = t82 * t76;
	t68 = t79 * t85 + t89;
	t66 = 0.1e1 / t68 ^ 2;
	t86 = t83 * t76;
	t87 = t82 * t78;
	t67 = t79 * t86 - t87;
	t94 = t66 * t67;
	t93 = t69 * t79;
	t73 = t77 ^ 2;
	t92 = t73 / t79 ^ 2;
	t91 = t77 * t83;
	t72 = 0.1e1 / (t82 ^ 2 * t92 + 0.1e1);
	t90 = t82 * t72;
	t84 = t67 ^ 2 * t66 + 0.1e1;
	t74 = 0.1e1 / t79;
	t65 = 0.1e1 / t68;
	t64 = (0.1e1 + t92) * t90;
	t63 = 0.1e1 / t84;
	t60 = 0.1e1 / t62;
	t59 = 0.1e1 / (t73 * t95 + 0.1e1);
	t58 = t84 * t63;
	t57 = (-t65 * t76 + t78 * t94) * t63 * t91;
	t56 = (t79 * t60 - (-t82 * t93 + t70 * t77 + (-t70 * t88 + t93) * t64) * t77 * t61) * t83 * t59;
	t1 = [t74 * t72 * t91, t64, t64, 0, 0, 0; (-t60 * t88 - (-t70 * t73 * t74 * t90 + (t72 - 0.1e1) * t77 * t69) * t77 * t95) * t59, t56, t56, 0, 0, 0; ((-t79 * t89 - t85) * t65 - (-t79 * t87 + t86) * t94) * t63, t57, t57, t58, t58, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end