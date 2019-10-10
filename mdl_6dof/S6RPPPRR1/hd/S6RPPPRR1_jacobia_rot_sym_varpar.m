% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPPRR1
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
%   Wie in S6RPPPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPPPRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_jacobia_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPPRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:25:52
	% EndTime: 2019-10-09 23:25:52
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (172->19), mult. (224->57), div. (52->9), fcn. (332->9), ass. (0->35)
	t37 = qJ(1) + pkin(9);
	t36 = cos(t37);
	t55 = t36 ^ 2;
	t42 = sin(qJ(5));
	t35 = sin(t37);
	t44 = cos(qJ(5));
	t50 = t35 * t44;
	t34 = atan2(t50, t42);
	t31 = sin(t34);
	t32 = cos(t34);
	t26 = t31 * t50 + t32 * t42;
	t25 = 0.1e1 / t26 ^ 2;
	t54 = t25 * t44;
	t41 = sin(qJ(6));
	t43 = cos(qJ(6));
	t46 = t42 * t43;
	t30 = -t35 * t41 + t36 * t46;
	t28 = 0.1e1 / t30 ^ 2;
	t47 = t41 * t42;
	t29 = t35 * t43 + t36 * t47;
	t53 = t28 * t29;
	t52 = t31 * t42;
	t40 = t44 ^ 2;
	t48 = 0.1e1 / t42 ^ 2 * t40;
	t33 = 0.1e1 / (t35 ^ 2 * t48 + 0.1e1);
	t51 = t35 * t33;
	t49 = t36 * t44;
	t45 = t29 ^ 2 * t28 + 0.1e1;
	t38 = 0.1e1 / t42;
	t27 = 0.1e1 / t30;
	t24 = 0.1e1 / t26;
	t23 = (-0.1e1 - t48) * t51;
	t22 = 0.1e1 / t45;
	t21 = 0.1e1 / (t55 * t40 * t25 + 0.1e1);
	t1 = [t38 * t33 * t49, 0, 0, 0, t23, 0; (t24 * t50 + (t32 * t38 * t40 * t51 + (-t33 + 0.1e1) * t44 * t31) * t55 * t54) * t21, 0, 0, 0, (t42 * t24 + (-t35 * t52 + t32 * t44 + (t32 * t50 - t52) * t23) * t54) * t36 * t21, 0; ((-t35 * t47 + t36 * t43) * t27 - (-t35 * t46 - t36 * t41) * t53) * t22, 0, 0, 0, (t27 * t41 - t43 * t53) * t22 * t49, t45 * t22;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,6);
end