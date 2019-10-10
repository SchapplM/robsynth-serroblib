% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR8
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
%   Wie in S6RPRRRR8_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRRR8_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR8_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (323->20), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t63 = qJ(3) + qJ(4);
	t61 = sin(t63);
	t62 = cos(t63);
	t67 = cos(qJ(1));
	t71 = t67 * t62;
	t56 = atan2(-t71, t61);
	t54 = sin(t56);
	t55 = cos(t56);
	t47 = -t54 * t71 + t55 * t61;
	t46 = 0.1e1 / t47 ^ 2;
	t65 = sin(qJ(1));
	t79 = t46 * t65 ^ 2;
	t64 = sin(qJ(5));
	t70 = t67 * t64;
	t66 = cos(qJ(5));
	t73 = t65 * t66;
	t53 = t61 * t73 + t70;
	t51 = 0.1e1 / t53 ^ 2;
	t69 = t67 * t66;
	t74 = t65 * t64;
	t52 = t61 * t74 - t69;
	t78 = t51 * t52;
	t77 = t54 * t61;
	t60 = t62 ^ 2;
	t76 = 0.1e1 / t61 ^ 2 * t60;
	t75 = t62 * t65;
	t57 = 0.1e1 / (t67 ^ 2 * t76 + 0.1e1);
	t72 = t67 * t57;
	t68 = t52 ^ 2 * t51 + 0.1e1;
	t58 = 0.1e1 / t61;
	t50 = 0.1e1 / t53;
	t49 = 0.1e1 / t68;
	t48 = (0.1e1 + t76) * t72;
	t45 = 0.1e1 / t47;
	t44 = 0.1e1 / (t60 * t79 + 0.1e1);
	t43 = (t50 * t64 - t66 * t78) * t49 * t75;
	t42 = (t61 * t45 + (t67 * t77 + t55 * t62 + (-t55 * t71 - t77) * t48) * t62 * t46) * t65 * t44;
	t1 = [t58 * t57 * t75, 0, t48, t48, 0, 0; (-t45 * t71 + (-t55 * t58 * t60 * t72 + (-t57 + 0.1e1) * t62 * t54) * t62 * t79) * t44, 0, t42, t42, 0, 0; ((t61 * t70 + t73) * t50 - (t61 * t69 - t74) * t78) * t49, 0, t43, t43, t68 * t49, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:07:54
	% EndTime: 2019-10-10 09:07:54
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (418->21), mult. (359->54), div. (84->9), fcn. (524->9), ass. (0->40)
	t81 = qJ(3) + qJ(4);
	t77 = sin(t81);
	t79 = cos(t81);
	t83 = cos(qJ(1));
	t85 = t83 * t79;
	t71 = atan2(-t85, t77);
	t69 = sin(t71);
	t70 = cos(t71);
	t62 = -t69 * t85 + t70 * t77;
	t61 = 0.1e1 / t62 ^ 2;
	t82 = sin(qJ(1));
	t95 = t61 * t82 ^ 2;
	t80 = qJ(5) + qJ(6);
	t76 = sin(t80);
	t87 = t83 * t76;
	t78 = cos(t80);
	t89 = t82 * t78;
	t68 = t77 * t89 + t87;
	t66 = 0.1e1 / t68 ^ 2;
	t86 = t83 * t78;
	t90 = t82 * t76;
	t67 = t77 * t90 - t86;
	t94 = t66 * t67;
	t93 = t69 * t77;
	t75 = t79 ^ 2;
	t92 = 0.1e1 / t77 ^ 2 * t75;
	t91 = t79 * t82;
	t72 = 0.1e1 / (t83 ^ 2 * t92 + 0.1e1);
	t88 = t83 * t72;
	t84 = t67 ^ 2 * t66 + 0.1e1;
	t73 = 0.1e1 / t77;
	t65 = 0.1e1 / t68;
	t64 = (0.1e1 + t92) * t88;
	t63 = 0.1e1 / t84;
	t60 = 0.1e1 / t62;
	t59 = 0.1e1 / (t75 * t95 + 0.1e1);
	t58 = t84 * t63;
	t57 = (t65 * t76 - t78 * t94) * t63 * t91;
	t56 = (t77 * t60 + (t83 * t93 + t70 * t79 + (-t70 * t85 - t93) * t64) * t79 * t61) * t82 * t59;
	t1 = [t73 * t72 * t91, 0, t64, t64, 0, 0; (-t60 * t85 + (-t70 * t73 * t75 * t88 + (-t72 + 0.1e1) * t79 * t69) * t79 * t95) * t59, 0, t56, t56, 0, 0; ((t77 * t87 + t89) * t65 - (t77 * t86 - t90) * t94) * t63, 0, t57, t57, t58, t58;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end