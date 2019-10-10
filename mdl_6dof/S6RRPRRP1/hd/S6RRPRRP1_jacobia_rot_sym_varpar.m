% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP1
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
%   Wie in S6RRPRRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:30
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRRP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:05
	% EndTime: 2019-10-10 10:30:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:05
	% EndTime: 2019-10-10 10:30:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:05
	% EndTime: 2019-10-10 10:30:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:05
	% EndTime: 2019-10-10 10:30:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:05
	% EndTime: 2019-10-10 10:30:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:06
	% EndTime: 2019-10-10 10:30:06
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (554->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t71 = qJ(2) + pkin(10) + qJ(4);
	t70 = cos(t71);
	t69 = sin(t71);
	t73 = sin(qJ(1));
	t81 = t73 * t69;
	t62 = atan2(-t81, -t70);
	t58 = sin(t62);
	t59 = cos(t62);
	t55 = -t58 * t81 - t59 * t70;
	t54 = 0.1e1 / t55 ^ 2;
	t75 = cos(qJ(1));
	t87 = t54 * t75 ^ 2;
	t86 = t58 * t70;
	t74 = cos(qJ(5));
	t77 = t75 * t74;
	t72 = sin(qJ(5));
	t80 = t73 * t72;
	t64 = t70 * t77 + t80;
	t61 = 0.1e1 / t64 ^ 2;
	t78 = t75 * t72;
	t79 = t73 * t74;
	t63 = t70 * t78 - t79;
	t85 = t61 * t63;
	t66 = t69 ^ 2;
	t84 = t66 / t70 ^ 2;
	t83 = t69 * t75;
	t65 = 0.1e1 / (t73 ^ 2 * t84 + 0.1e1);
	t82 = t73 * t65;
	t76 = t63 ^ 2 * t61 + 0.1e1;
	t67 = 0.1e1 / t70;
	t60 = 0.1e1 / t64;
	t57 = 0.1e1 / t76;
	t56 = (0.1e1 + t84) * t82;
	t53 = 0.1e1 / t55;
	t52 = 0.1e1 / (t66 * t87 + 0.1e1);
	t51 = (-t60 * t72 + t74 * t85) * t57 * t83;
	t50 = (t70 * t53 - (-t73 * t86 + t59 * t69 + (-t59 * t81 + t86) * t56) * t69 * t54) * t75 * t52;
	t1 = [t67 * t65 * t83, t56, 0, t56, 0, 0; (-t53 * t81 - (-t59 * t66 * t67 * t82 + (t65 - 0.1e1) * t69 * t58) * t69 * t87) * t52, t50, 0, t50, 0, 0; ((-t70 * t80 - t77) * t60 - (-t70 * t79 + t78) * t85) * t57, t51, 0, t51, t76 * t57, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:30:06
	% EndTime: 2019-10-10 10:30:06
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (554->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t69 = qJ(2) + pkin(10) + qJ(4);
	t68 = cos(t69);
	t67 = sin(t69);
	t71 = sin(qJ(1));
	t79 = t71 * t67;
	t60 = atan2(-t79, -t68);
	t56 = sin(t60);
	t57 = cos(t60);
	t53 = -t56 * t79 - t57 * t68;
	t52 = 0.1e1 / t53 ^ 2;
	t73 = cos(qJ(1));
	t85 = t52 * t73 ^ 2;
	t84 = t56 * t68;
	t72 = cos(qJ(5));
	t75 = t73 * t72;
	t70 = sin(qJ(5));
	t78 = t71 * t70;
	t62 = t68 * t75 + t78;
	t59 = 0.1e1 / t62 ^ 2;
	t76 = t73 * t70;
	t77 = t71 * t72;
	t61 = t68 * t76 - t77;
	t83 = t59 * t61;
	t64 = t67 ^ 2;
	t82 = t64 / t68 ^ 2;
	t81 = t67 * t73;
	t63 = 0.1e1 / (t71 ^ 2 * t82 + 0.1e1);
	t80 = t71 * t63;
	t74 = t61 ^ 2 * t59 + 0.1e1;
	t65 = 0.1e1 / t68;
	t58 = 0.1e1 / t62;
	t55 = 0.1e1 / t74;
	t54 = (0.1e1 + t82) * t80;
	t51 = 0.1e1 / t53;
	t50 = 0.1e1 / (t64 * t85 + 0.1e1);
	t49 = (-t58 * t70 + t72 * t83) * t55 * t81;
	t48 = (t68 * t51 - (-t71 * t84 + t57 * t67 + (-t57 * t79 + t84) * t54) * t67 * t52) * t73 * t50;
	t1 = [t65 * t63 * t81, t54, 0, t54, 0, 0; (-t51 * t79 - (-t57 * t64 * t65 * t80 + (t63 - 0.1e1) * t67 * t56) * t67 * t85) * t50, t48, 0, t48, 0, 0; ((-t68 * t78 - t75) * t58 - (-t68 * t77 + t76) * t83) * t55, t49, 0, t49, t74 * t55, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end