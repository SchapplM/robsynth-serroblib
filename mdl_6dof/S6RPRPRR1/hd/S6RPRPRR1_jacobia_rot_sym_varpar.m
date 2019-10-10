% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR1
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
%   Wie in S6RPRPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_jacobia_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:46:14
	% EndTime: 2019-10-10 00:46:14
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (713->22), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->39)
	t74 = qJ(3) + pkin(11) + qJ(5);
	t71 = cos(t74);
	t70 = sin(t74);
	t75 = qJ(1) + pkin(10);
	t72 = sin(t75);
	t83 = t72 * t70;
	t65 = atan2(-t83, -t71);
	t63 = sin(t65);
	t64 = cos(t65);
	t56 = -t63 * t83 - t64 * t71;
	t55 = 0.1e1 / t56 ^ 2;
	t73 = cos(t75);
	t89 = t55 * t73 ^ 2;
	t77 = cos(qJ(6));
	t79 = t73 * t77;
	t76 = sin(qJ(6));
	t82 = t72 * t76;
	t62 = t71 * t79 + t82;
	t60 = 0.1e1 / t62 ^ 2;
	t80 = t73 * t76;
	t81 = t72 * t77;
	t61 = t71 * t80 - t81;
	t88 = t60 * t61;
	t87 = t63 * t71;
	t67 = t70 ^ 2;
	t86 = t67 / t71 ^ 2;
	t85 = t70 * t73;
	t66 = 0.1e1 / (t72 ^ 2 * t86 + 0.1e1);
	t84 = t72 * t66;
	t78 = t61 ^ 2 * t60 + 0.1e1;
	t68 = 0.1e1 / t71;
	t59 = 0.1e1 / t62;
	t58 = 0.1e1 / t78;
	t57 = (0.1e1 + t86) * t84;
	t54 = 0.1e1 / t56;
	t53 = 0.1e1 / (t67 * t89 + 0.1e1);
	t52 = (-t59 * t76 + t77 * t88) * t58 * t85;
	t51 = (t71 * t54 - (-t72 * t87 + t64 * t70 + (-t64 * t83 + t87) * t57) * t70 * t55) * t73 * t53;
	t1 = [t68 * t66 * t85, 0, t57, 0, t57, 0; (-t54 * t83 - (-t64 * t67 * t68 * t84 + (t66 - 0.1e1) * t70 * t63) * t70 * t89) * t53, 0, t51, 0, t51, 0; ((-t71 * t82 - t79) * t59 - (-t71 * t81 + t80) * t88) * t58, 0, t52, 0, t52, t78 * t58;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end