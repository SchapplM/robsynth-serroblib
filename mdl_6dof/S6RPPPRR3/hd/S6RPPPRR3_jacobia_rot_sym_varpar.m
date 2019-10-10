% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR3
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
%   Wie in S6RPPPRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPPRR3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobia_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->0), mult. (22->0), div. (6->0), fcn. (32->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:29:15
	% EndTime: 2019-10-09 23:29:15
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (325->21), mult. (442->61), div. (52->9), fcn. (659->11), ass. (0->39)
	t59 = sin(pkin(9));
	t60 = cos(pkin(9));
	t70 = sin(qJ(1));
	t71 = cos(qJ(1));
	t43 = -t70 * t59 - t71 * t60;
	t72 = t43 ^ 2;
	t55 = pkin(10) + qJ(5);
	t54 = cos(t55);
	t44 = t71 * t59 - t70 * t60;
	t53 = sin(t55);
	t64 = t44 * t53;
	t41 = atan2(t64, t54);
	t39 = sin(t41);
	t40 = cos(t41);
	t33 = t39 * t64 + t40 * t54;
	t32 = 0.1e1 / t33 ^ 2;
	t69 = t32 * t53;
	t56 = sin(qJ(6));
	t57 = cos(qJ(6));
	t61 = t54 * t57;
	t38 = -t43 * t61 + t44 * t56;
	t36 = 0.1e1 / t38 ^ 2;
	t62 = t54 * t56;
	t37 = -t43 * t62 - t44 * t57;
	t68 = t36 * t37;
	t67 = t39 * t54;
	t66 = t43 * t53;
	t50 = t53 ^ 2;
	t63 = t50 / t54 ^ 2;
	t42 = 0.1e1 / (t44 ^ 2 * t63 + 0.1e1);
	t65 = t44 * t42;
	t58 = t37 ^ 2 * t36 + 0.1e1;
	t51 = 0.1e1 / t54;
	t35 = 0.1e1 / t38;
	t34 = 0.1e1 / t58;
	t31 = 0.1e1 / t33;
	t30 = (0.1e1 + t63) * t65;
	t29 = 0.1e1 / (t72 * t50 * t32 + 0.1e1);
	t1 = [t51 * t42 * t66, 0, 0, 0, t30, 0; (t31 * t64 + (t40 * t50 * t51 * t65 + (-t42 + 0.1e1) * t53 * t39) * t72 * t69) * t29, 0, 0, 0, (-t54 * t31 + (t44 * t67 - t40 * t53 + (t40 * t64 - t67) * t30) * t69) * t43 * t29, 0; ((-t43 * t57 + t44 * t62) * t35 - (t43 * t56 + t44 * t61) * t68) * t34, 0, 0, 0, (t35 * t56 - t57 * t68) * t34 * t66, t58 * t34;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end