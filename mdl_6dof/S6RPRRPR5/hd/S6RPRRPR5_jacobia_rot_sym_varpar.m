% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR5
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
%   Wie in S6RPRRPR5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPR5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (467->18), mult. (231->39), div. (67->9), fcn. (349->7), ass. (0->30)
	t56 = cos(qJ(1));
	t54 = t56 ^ 2;
	t51 = pkin(10) + qJ(3) + qJ(4);
	t50 = cos(t51);
	t49 = sin(t51);
	t55 = sin(qJ(1));
	t60 = t55 * t49;
	t43 = atan2(-t60, -t50);
	t41 = sin(t43);
	t42 = cos(t43);
	t39 = -t41 * t60 - t42 * t50;
	t38 = 0.1e1 / t39 ^ 2;
	t66 = t38 * t49;
	t65 = t41 * t50;
	t46 = t49 ^ 2;
	t57 = t50 ^ 2;
	t64 = t46 / t57;
	t63 = t49 * t56;
	t58 = t55 ^ 2;
	t62 = 0.1e1 / t58 * t54;
	t44 = 0.1e1 / (t58 * t64 + 0.1e1);
	t61 = t55 * t44;
	t45 = 0.1e1 / (t57 * t62 + 0.1e1);
	t59 = 0.1e1 / t55 * t45 * t63;
	t47 = 0.1e1 / t50;
	t40 = (0.1e1 + t64) * t61;
	t37 = 0.1e1 / t39;
	t36 = 0.1e1 / (t54 * t46 * t38 + 0.1e1);
	t35 = (t50 * t37 - (-t55 * t65 + t42 * t49 + (-t42 * t60 + t65) * t40) * t66) * t56 * t36;
	t1 = [t47 * t44 * t63, 0, t40, t40, 0, 0; (-t37 * t60 - (-t42 * t46 * t47 * t61 + (t44 - 0.1e1) * t49 * t41) * t54 * t66) * t36, 0, t35, t35, 0, 0; (-0.1e1 - t62) * t50 * t45, 0, -t59, -t59, 0, 0;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:30:29
	% EndTime: 2019-10-10 01:30:29
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (520->21), mult. (332->54), div. (79->9), fcn. (489->9), ass. (0->38)
	t71 = pkin(10) + qJ(3) + qJ(4);
	t69 = sin(t71);
	t70 = cos(t71);
	t73 = sin(qJ(1));
	t81 = t73 * t70;
	t64 = atan2(-t81, t69);
	t58 = sin(t64);
	t59 = cos(t64);
	t55 = -t58 * t81 + t59 * t69;
	t54 = 0.1e1 / t55 ^ 2;
	t75 = cos(qJ(1));
	t87 = t54 * t75 ^ 2;
	t86 = t58 * t69;
	t72 = sin(qJ(6));
	t78 = t75 * t72;
	t74 = cos(qJ(6));
	t79 = t73 * t74;
	t63 = t69 * t78 + t79;
	t61 = 0.1e1 / t63 ^ 2;
	t77 = t75 * t74;
	t80 = t73 * t72;
	t62 = -t69 * t77 + t80;
	t85 = t61 * t62;
	t68 = t70 ^ 2;
	t84 = 0.1e1 / t69 ^ 2 * t68;
	t83 = t70 * t75;
	t65 = 0.1e1 / (t73 ^ 2 * t84 + 0.1e1);
	t82 = t73 * t65;
	t76 = t62 ^ 2 * t61 + 0.1e1;
	t66 = 0.1e1 / t69;
	t60 = 0.1e1 / t63;
	t57 = 0.1e1 / t76;
	t56 = (0.1e1 + t84) * t82;
	t53 = 0.1e1 / t55;
	t52 = 0.1e1 / (t68 * t87 + 0.1e1);
	t51 = (-t60 * t74 - t72 * t85) * t57 * t83;
	t50 = (-t69 * t53 - (t73 * t86 + t59 * t70 + (-t59 * t81 - t86) * t56) * t70 * t54) * t75 * t52;
	t1 = [-t66 * t65 * t83, 0, t56, t56, 0, 0; (-t53 * t81 - (t59 * t66 * t68 * t82 + (t65 - 0.1e1) * t70 * t58) * t70 * t87) * t52, 0, t50, t50, 0, 0; ((t69 * t79 + t78) * t60 - (-t69 * t80 + t77) * t85) * t57, 0, t51, t51, 0, t76 * t57;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end