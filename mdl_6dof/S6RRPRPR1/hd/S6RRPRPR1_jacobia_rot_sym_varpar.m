% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR1
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
%   Wie in S6RRPRPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:04
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:25
	% EndTime: 2019-10-10 10:04:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:25
	% EndTime: 2019-10-10 10:04:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:25
	% EndTime: 2019-10-10 10:04:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:25
	% EndTime: 2019-10-10 10:04:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:25
	% EndTime: 2019-10-10 10:04:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (530->21), mult. (305->53), div. (74->9), fcn. (454->9), ass. (0->37)
	t65 = qJ(2) + pkin(10) + qJ(4);
	t64 = cos(t65);
	t63 = sin(t65);
	t68 = sin(qJ(1));
	t74 = t68 * t63;
	t58 = atan2(-t74, -t64);
	t52 = sin(t58);
	t53 = cos(t58);
	t49 = -t52 * t74 - t53 * t64;
	t48 = 0.1e1 / t49 ^ 2;
	t69 = cos(qJ(1));
	t80 = t48 * t69 ^ 2;
	t79 = t52 * t64;
	t67 = cos(pkin(11));
	t70 = t69 * t67;
	t66 = sin(pkin(11));
	t73 = t68 * t66;
	t57 = t64 * t70 + t73;
	t55 = 0.1e1 / t57 ^ 2;
	t71 = t69 * t66;
	t72 = t68 * t67;
	t56 = t64 * t71 - t72;
	t78 = t55 * t56;
	t60 = t63 ^ 2;
	t77 = t60 / t64 ^ 2;
	t76 = t63 * t69;
	t59 = 0.1e1 / (t68 ^ 2 * t77 + 0.1e1);
	t75 = t68 * t59;
	t61 = 0.1e1 / t64;
	t54 = 0.1e1 / t57;
	t51 = 0.1e1 / (t56 ^ 2 * t55 + 0.1e1);
	t50 = (0.1e1 + t77) * t75;
	t47 = 0.1e1 / t49;
	t46 = 0.1e1 / (t60 * t80 + 0.1e1);
	t45 = (-t54 * t66 + t67 * t78) * t51 * t76;
	t44 = (t64 * t47 - (-t68 * t79 + t53 * t63 + (-t53 * t74 + t79) * t50) * t63 * t48) * t69 * t46;
	t1 = [t61 * t59 * t76, t50, 0, t50, 0, 0; (-t47 * t74 - (-t53 * t60 * t61 * t75 + (t59 - 0.1e1) * t63 * t52) * t63 * t80) * t46, t44, 0, t44, 0, 0; ((-t64 * t73 - t70) * t54 - (-t64 * t72 + t71) * t78) * t51, t45, 0, t45, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:04:26
	% EndTime: 2019-10-10 10:04:26
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (618->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
	t72 = qJ(2) + pkin(10) + qJ(4);
	t69 = cos(t72);
	t68 = sin(t72);
	t74 = sin(qJ(1));
	t81 = t74 * t68;
	t63 = atan2(-t81, -t69);
	t61 = sin(t63);
	t62 = cos(t63);
	t54 = -t61 * t81 - t62 * t69;
	t53 = 0.1e1 / t54 ^ 2;
	t75 = cos(qJ(1));
	t87 = t53 * t75 ^ 2;
	t73 = pkin(11) + qJ(6);
	t71 = cos(t73);
	t77 = t75 * t71;
	t70 = sin(t73);
	t80 = t74 * t70;
	t60 = t69 * t77 + t80;
	t58 = 0.1e1 / t60 ^ 2;
	t78 = t75 * t70;
	t79 = t74 * t71;
	t59 = t69 * t78 - t79;
	t86 = t58 * t59;
	t85 = t61 * t69;
	t65 = t68 ^ 2;
	t84 = t65 / t69 ^ 2;
	t83 = t68 * t75;
	t64 = 0.1e1 / (t74 ^ 2 * t84 + 0.1e1);
	t82 = t74 * t64;
	t76 = t59 ^ 2 * t58 + 0.1e1;
	t66 = 0.1e1 / t69;
	t57 = 0.1e1 / t60;
	t56 = 0.1e1 / t76;
	t55 = (0.1e1 + t84) * t82;
	t52 = 0.1e1 / t54;
	t51 = 0.1e1 / (t65 * t87 + 0.1e1);
	t50 = (-t57 * t70 + t71 * t86) * t56 * t83;
	t49 = (t69 * t52 - (-t74 * t85 + t62 * t68 + (-t62 * t81 + t85) * t55) * t68 * t53) * t75 * t51;
	t1 = [t66 * t64 * t83, t55, 0, t55, 0, 0; (-t52 * t81 - (-t62 * t65 * t66 * t82 + (t64 - 0.1e1) * t68 * t61) * t68 * t87) * t51, t49, 0, t49, 0, 0; ((-t69 * t80 - t77) * t57 - (-t69 * t79 + t78) * t86) * t56, t50, 0, t50, 0, t76 * t56;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end