% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP1
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
%   Wie in S6RRRPRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPRP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:35:10
	% EndTime: 2019-10-10 11:35:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:35:10
	% EndTime: 2019-10-10 11:35:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:35:10
	% EndTime: 2019-10-10 11:35:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:35:10
	% EndTime: 2019-10-10 11:35:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:35:10
	% EndTime: 2019-10-10 11:35:10
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:35:10
	% EndTime: 2019-10-10 11:35:10
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (554->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t74 = qJ(2) + qJ(3) + pkin(10);
	t73 = cos(t74);
	t72 = sin(t74);
	t76 = sin(qJ(1));
	t84 = t76 * t72;
	t65 = atan2(-t84, -t73);
	t61 = sin(t65);
	t62 = cos(t65);
	t58 = -t61 * t84 - t62 * t73;
	t57 = 0.1e1 / t58 ^ 2;
	t78 = cos(qJ(1));
	t90 = t57 * t78 ^ 2;
	t89 = t61 * t73;
	t77 = cos(qJ(5));
	t80 = t78 * t77;
	t75 = sin(qJ(5));
	t83 = t76 * t75;
	t67 = t73 * t80 + t83;
	t64 = 0.1e1 / t67 ^ 2;
	t81 = t78 * t75;
	t82 = t76 * t77;
	t66 = t73 * t81 - t82;
	t88 = t64 * t66;
	t69 = t72 ^ 2;
	t87 = t69 / t73 ^ 2;
	t86 = t72 * t78;
	t68 = 0.1e1 / (t76 ^ 2 * t87 + 0.1e1);
	t85 = t76 * t68;
	t79 = t66 ^ 2 * t64 + 0.1e1;
	t70 = 0.1e1 / t73;
	t63 = 0.1e1 / t67;
	t60 = 0.1e1 / t79;
	t59 = (0.1e1 + t87) * t85;
	t56 = 0.1e1 / t58;
	t55 = 0.1e1 / (t69 * t90 + 0.1e1);
	t54 = (-t63 * t75 + t77 * t88) * t60 * t86;
	t53 = (t73 * t56 - (-t76 * t89 + t62 * t72 + (-t62 * t84 + t89) * t59) * t72 * t57) * t78 * t55;
	t1 = [t70 * t68 * t86, t59, t59, 0, 0, 0; (-t56 * t84 - (-t62 * t69 * t70 * t85 + (t68 - 0.1e1) * t72 * t61) * t72 * t90) * t55, t53, t53, 0, 0, 0; ((-t73 * t83 - t80) * t63 - (-t73 * t82 + t81) * t88) * t60, t54, t54, 0, t79 * t60, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:35:10
	% EndTime: 2019-10-10 11:35:10
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (554->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t72 = qJ(2) + qJ(3) + pkin(10);
	t71 = cos(t72);
	t70 = sin(t72);
	t74 = sin(qJ(1));
	t82 = t74 * t70;
	t63 = atan2(-t82, -t71);
	t59 = sin(t63);
	t60 = cos(t63);
	t56 = -t59 * t82 - t60 * t71;
	t55 = 0.1e1 / t56 ^ 2;
	t76 = cos(qJ(1));
	t88 = t55 * t76 ^ 2;
	t87 = t59 * t71;
	t75 = cos(qJ(5));
	t78 = t76 * t75;
	t73 = sin(qJ(5));
	t81 = t74 * t73;
	t65 = t71 * t78 + t81;
	t62 = 0.1e1 / t65 ^ 2;
	t79 = t76 * t73;
	t80 = t74 * t75;
	t64 = t71 * t79 - t80;
	t86 = t62 * t64;
	t67 = t70 ^ 2;
	t85 = t67 / t71 ^ 2;
	t84 = t70 * t76;
	t66 = 0.1e1 / (t74 ^ 2 * t85 + 0.1e1);
	t83 = t74 * t66;
	t77 = t64 ^ 2 * t62 + 0.1e1;
	t68 = 0.1e1 / t71;
	t61 = 0.1e1 / t65;
	t58 = 0.1e1 / t77;
	t57 = (0.1e1 + t85) * t83;
	t54 = 0.1e1 / t56;
	t53 = 0.1e1 / (t67 * t88 + 0.1e1);
	t52 = (-t61 * t73 + t75 * t86) * t58 * t84;
	t51 = (t71 * t54 - (-t74 * t87 + t60 * t70 + (-t60 * t82 + t87) * t57) * t70 * t55) * t76 * t53;
	t1 = [t68 * t66 * t84, t57, t57, 0, 0, 0; (-t54 * t82 - (-t60 * t67 * t68 * t83 + (t66 - 0.1e1) * t70 * t59) * t70 * t88) * t53, t51, t51, 0, 0, 0; ((-t71 * t81 - t78) * t61 - (-t71 * t80 + t79) * t86) * t58, t52, t52, 0, t77 * t58, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end