% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRR1
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
%   Wie in S6RPPRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:01
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPRRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:01:23
	% EndTime: 2019-10-10 00:01:23
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (713->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
	t72 = pkin(11) + qJ(4) + qJ(5);
	t69 = cos(t72);
	t68 = sin(t72);
	t73 = qJ(1) + pkin(10);
	t70 = sin(t73);
	t81 = t70 * t68;
	t63 = atan2(-t81, -t69);
	t61 = sin(t63);
	t62 = cos(t63);
	t54 = -t61 * t81 - t62 * t69;
	t53 = 0.1e1 / t54 ^ 2;
	t71 = cos(t73);
	t87 = t53 * t71 ^ 2;
	t75 = cos(qJ(6));
	t77 = t71 * t75;
	t74 = sin(qJ(6));
	t80 = t70 * t74;
	t60 = t69 * t77 + t80;
	t58 = 0.1e1 / t60 ^ 2;
	t78 = t71 * t74;
	t79 = t70 * t75;
	t59 = t69 * t78 - t79;
	t86 = t58 * t59;
	t85 = t61 * t69;
	t65 = t68 ^ 2;
	t84 = t65 / t69 ^ 2;
	t83 = t68 * t71;
	t64 = 0.1e1 / (t70 ^ 2 * t84 + 0.1e1);
	t82 = t70 * t64;
	t76 = t59 ^ 2 * t58 + 0.1e1;
	t66 = 0.1e1 / t69;
	t57 = 0.1e1 / t60;
	t56 = 0.1e1 / t76;
	t55 = (0.1e1 + t84) * t82;
	t52 = 0.1e1 / t54;
	t51 = 0.1e1 / (t65 * t87 + 0.1e1);
	t50 = (-t57 * t74 + t75 * t86) * t56 * t83;
	t49 = (t69 * t52 - (-t70 * t85 + t62 * t68 + (-t62 * t81 + t85) * t55) * t68 * t53) * t71 * t51;
	t1 = [t66 * t64 * t83, 0, 0, t55, t55, 0; (-t52 * t81 - (-t62 * t65 * t66 * t82 + (t64 - 0.1e1) * t68 * t61) * t68 * t87) * t51, 0, 0, t49, t49, 0; ((-t69 * t80 - t77) * t57 - (-t69 * t79 + t78) * t86) * t56, 0, 0, t50, t50, t76 * t56;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end