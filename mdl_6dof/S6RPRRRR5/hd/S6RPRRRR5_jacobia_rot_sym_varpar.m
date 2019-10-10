% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR5
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
%   Wie in S6RPRRRR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:06
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (554->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t69 = pkin(11) + qJ(3) + qJ(4);
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
	t1 = [t65 * t63 * t81, 0, t54, t54, 0, 0; (-t51 * t79 - (-t57 * t64 * t65 * t80 + (t63 - 0.1e1) * t67 * t56) * t67 * t85) * t50, 0, t48, t48, 0, 0; ((-t68 * t78 - t75) * t58 - (-t68 * t77 + t76) * t83) * t55, 0, t49, t49, t74 * t55, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:06:07
	% EndTime: 2019-10-10 09:06:07
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (656->22), mult. (359->54), div. (84->9), fcn. (524->9), ass. (0->40)
	t78 = pkin(11) + qJ(3) + qJ(4);
	t77 = cos(t78);
	t76 = sin(t78);
	t82 = sin(qJ(1));
	t89 = t82 * t76;
	t71 = atan2(-t89, -t77);
	t69 = sin(t71);
	t70 = cos(t71);
	t62 = -t69 * t89 - t70 * t77;
	t61 = 0.1e1 / t62 ^ 2;
	t83 = cos(qJ(1));
	t95 = t61 * t83 ^ 2;
	t81 = qJ(5) + qJ(6);
	t80 = cos(t81);
	t85 = t83 * t80;
	t79 = sin(t81);
	t88 = t82 * t79;
	t68 = t77 * t85 + t88;
	t66 = 0.1e1 / t68 ^ 2;
	t86 = t83 * t79;
	t87 = t82 * t80;
	t67 = t77 * t86 - t87;
	t94 = t66 * t67;
	t93 = t69 * t77;
	t73 = t76 ^ 2;
	t92 = t73 / t77 ^ 2;
	t91 = t76 * t83;
	t72 = 0.1e1 / (t82 ^ 2 * t92 + 0.1e1);
	t90 = t82 * t72;
	t84 = t67 ^ 2 * t66 + 0.1e1;
	t74 = 0.1e1 / t77;
	t65 = 0.1e1 / t68;
	t64 = 0.1e1 / t84;
	t63 = (0.1e1 + t92) * t90;
	t60 = 0.1e1 / t62;
	t59 = 0.1e1 / (t73 * t95 + 0.1e1);
	t58 = t84 * t64;
	t57 = (-t65 * t79 + t80 * t94) * t64 * t91;
	t56 = (t77 * t60 - (-t82 * t93 + t70 * t76 + (-t70 * t89 + t93) * t63) * t76 * t61) * t83 * t59;
	t1 = [t74 * t72 * t91, 0, t63, t63, 0, 0; (-t60 * t89 - (-t70 * t73 * t74 * t90 + (t72 - 0.1e1) * t76 * t69) * t76 * t95) * t59, 0, t56, t56, 0, 0; ((-t77 * t88 - t85) * t65 - (-t77 * t87 + t86) * t94) * t64, 0, t57, t57, t58, t58;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end