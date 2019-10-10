% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR1
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
%   Wie in S6RRRPPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:17
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (530->21), mult. (305->53), div. (74->9), fcn. (454->9), ass. (0->37)
	t67 = qJ(2) + qJ(3) + pkin(10);
	t66 = cos(t67);
	t65 = sin(t67);
	t70 = sin(qJ(1));
	t76 = t70 * t65;
	t60 = atan2(-t76, -t66);
	t54 = sin(t60);
	t55 = cos(t60);
	t51 = -t54 * t76 - t55 * t66;
	t50 = 0.1e1 / t51 ^ 2;
	t71 = cos(qJ(1));
	t82 = t50 * t71 ^ 2;
	t81 = t54 * t66;
	t69 = cos(pkin(11));
	t72 = t71 * t69;
	t68 = sin(pkin(11));
	t75 = t70 * t68;
	t59 = t66 * t72 + t75;
	t57 = 0.1e1 / t59 ^ 2;
	t73 = t71 * t68;
	t74 = t70 * t69;
	t58 = t66 * t73 - t74;
	t80 = t57 * t58;
	t62 = t65 ^ 2;
	t79 = t62 / t66 ^ 2;
	t78 = t65 * t71;
	t61 = 0.1e1 / (t70 ^ 2 * t79 + 0.1e1);
	t77 = t70 * t61;
	t63 = 0.1e1 / t66;
	t56 = 0.1e1 / t59;
	t53 = 0.1e1 / (t58 ^ 2 * t57 + 0.1e1);
	t52 = (0.1e1 + t79) * t77;
	t49 = 0.1e1 / t51;
	t48 = 0.1e1 / (t62 * t82 + 0.1e1);
	t47 = (-t56 * t68 + t69 * t80) * t53 * t78;
	t46 = (t66 * t49 - (-t70 * t81 + t55 * t65 + (-t55 * t76 + t81) * t52) * t65 * t50) * t71 * t48;
	t1 = [t63 * t61 * t78, t52, t52, 0, 0, 0; (-t49 * t76 - (-t55 * t62 * t63 * t77 + (t61 - 0.1e1) * t65 * t54) * t65 * t82) * t48, t46, t46, 0, 0, 0; ((-t66 * t75 - t72) * t56 - (-t66 * t74 + t73) * t80) * t53, t47, t47, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:17:04
	% EndTime: 2019-10-10 11:17:04
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (618->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
	t75 = qJ(2) + qJ(3) + pkin(10);
	t72 = cos(t75);
	t71 = sin(t75);
	t77 = sin(qJ(1));
	t84 = t77 * t71;
	t66 = atan2(-t84, -t72);
	t64 = sin(t66);
	t65 = cos(t66);
	t57 = -t64 * t84 - t65 * t72;
	t56 = 0.1e1 / t57 ^ 2;
	t78 = cos(qJ(1));
	t90 = t56 * t78 ^ 2;
	t76 = pkin(11) + qJ(6);
	t74 = cos(t76);
	t80 = t78 * t74;
	t73 = sin(t76);
	t83 = t77 * t73;
	t63 = t72 * t80 + t83;
	t61 = 0.1e1 / t63 ^ 2;
	t81 = t78 * t73;
	t82 = t77 * t74;
	t62 = t72 * t81 - t82;
	t89 = t61 * t62;
	t88 = t64 * t72;
	t68 = t71 ^ 2;
	t87 = t68 / t72 ^ 2;
	t86 = t71 * t78;
	t67 = 0.1e1 / (t77 ^ 2 * t87 + 0.1e1);
	t85 = t77 * t67;
	t79 = t62 ^ 2 * t61 + 0.1e1;
	t69 = 0.1e1 / t72;
	t60 = 0.1e1 / t63;
	t59 = 0.1e1 / t79;
	t58 = (0.1e1 + t87) * t85;
	t55 = 0.1e1 / t57;
	t54 = 0.1e1 / (t68 * t90 + 0.1e1);
	t53 = (-t60 * t73 + t74 * t89) * t59 * t86;
	t52 = (t72 * t55 - (-t77 * t88 + t65 * t71 + (-t65 * t84 + t88) * t58) * t71 * t56) * t78 * t54;
	t1 = [t69 * t67 * t86, t58, t58, 0, 0, 0; (-t55 * t84 - (-t65 * t68 * t69 * t85 + (t67 - 0.1e1) * t71 * t64) * t71 * t90) * t54, t52, t52, 0, 0, 0; ((-t72 * t83 - t80) * t60 - (-t72 * t82 + t81) * t89) * t59, t53, t53, 0, 0, t79 * t59;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end