% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR7
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
%   Wie in S6RPRRPR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:55
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:33:56
	% EndTime: 2019-10-10 01:33:56
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (519->20), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t68 = qJ(3) + qJ(4) + pkin(10);
	t66 = sin(t68);
	t67 = cos(t68);
	t72 = cos(qJ(1));
	t76 = t72 * t67;
	t61 = atan2(-t76, t66);
	t55 = sin(t61);
	t56 = cos(t61);
	t52 = -t55 * t76 + t56 * t66;
	t51 = 0.1e1 / t52 ^ 2;
	t70 = sin(qJ(1));
	t84 = t51 * t70 ^ 2;
	t83 = t55 * t66;
	t69 = sin(qJ(6));
	t75 = t72 * t69;
	t71 = cos(qJ(6));
	t77 = t70 * t71;
	t60 = t66 * t77 + t75;
	t58 = 0.1e1 / t60 ^ 2;
	t74 = t72 * t71;
	t78 = t70 * t69;
	t59 = t66 * t78 - t74;
	t82 = t58 * t59;
	t65 = t67 ^ 2;
	t80 = 0.1e1 / t66 ^ 2 * t65;
	t62 = 0.1e1 / (t72 ^ 2 * t80 + 0.1e1);
	t81 = t62 * t72;
	t79 = t67 * t70;
	t73 = t59 ^ 2 * t58 + 0.1e1;
	t63 = 0.1e1 / t66;
	t57 = 0.1e1 / t60;
	t54 = 0.1e1 / t73;
	t53 = (0.1e1 + t80) * t81;
	t50 = 0.1e1 / t52;
	t49 = 0.1e1 / (t65 * t84 + 0.1e1);
	t48 = (t57 * t69 - t71 * t82) * t54 * t79;
	t47 = (t66 * t50 + (t72 * t83 + t56 * t67 + (-t56 * t76 - t83) * t53) * t67 * t51) * t70 * t49;
	t1 = [t63 * t62 * t79, 0, t53, t53, 0, 0; (-t50 * t76 + (-t56 * t63 * t65 * t81 + (-t62 + 0.1e1) * t67 * t55) * t67 * t84) * t49, 0, t47, t47, 0, 0; ((t66 * t75 + t77) * t57 - (t66 * t74 - t78) * t82) * t54, 0, t48, t48, 0, t73 * t54;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end