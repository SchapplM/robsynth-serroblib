% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRRR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:56
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRRR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:53
	% EndTime: 2019-12-29 20:55:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:54
	% EndTime: 2019-12-29 20:55:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:53
	% EndTime: 2019-12-29 20:55:54
	% DurationCPUTime: 0.25s
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
	t1 = [t64 * t62 * t80, t53, t53, 0, 0; (-t50 * t78 - (-t60 * t63 * t64 * t79 + (t62 - 0.1e1) * t66 * t59) * t66 * t84) * t49, t47, t47, 0, 0; ((-t67 * t77 - t74) * t55 - (-t67 * t76 + t75) * t83) * t54, t48, t48, t73 * t54, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:55:53
	% EndTime: 2019-12-29 20:55:53
	% DurationCPUTime: 0.26s
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
	t1 = [t73 * t71 * t90, t63, t63, 0, 0; (-t59 * t87 - (-t69 * t72 * t73 * t89 + (t71 - 0.1e1) * t76 * t68) * t76 * t94) * t58, t55, t55, 0, 0; ((-t78 * t88 - t84) * t64 - (-t78 * t86 + t85) * t93) * t62, t56, t56, t57, t57;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end