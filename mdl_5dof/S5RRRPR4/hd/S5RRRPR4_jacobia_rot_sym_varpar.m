% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 19:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRPR4_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR4_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:56:52
	% EndTime: 2019-12-29 19:56:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:56:53
	% EndTime: 2019-12-29 19:56:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:56:48
	% EndTime: 2019-12-29 19:56:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (18->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:56:53
	% EndTime: 2019-12-29 19:56:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:56:53
	% EndTime: 2019-12-29 19:56:53
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (223->17), mult. (208->39), div. (72->9), fcn. (333->7), ass. (0->30)
	t45 = qJ(1) + qJ(2);
	t41 = cos(t45);
	t48 = t41 ^ 2;
	t47 = cos(qJ(3));
	t40 = sin(t45);
	t46 = sin(qJ(3));
	t50 = t40 * t46;
	t36 = atan2(-t50, -t47);
	t33 = sin(t36);
	t34 = cos(t36);
	t30 = -t33 * t50 - t34 * t47;
	t29 = 0.1e1 / t30 ^ 2;
	t54 = t29 * t46;
	t53 = t33 * t47;
	t38 = t40 ^ 2;
	t52 = t38 / t48;
	t42 = t46 ^ 2;
	t44 = 0.1e1 / t47 ^ 2;
	t49 = t42 * t44;
	t37 = 0.1e1 / (t38 * t49 + 0.1e1);
	t51 = t40 * t37;
	t43 = 0.1e1 / t47;
	t35 = 0.1e1 / (t44 * t52 + 0.1e1);
	t32 = t41 * t46 * t43 * t37;
	t31 = (0.1e1 + t49) * t51;
	t28 = 0.1e1 / t30;
	t27 = (-0.1e1 - t52) * t43 * t35;
	t26 = 0.1e1 / (t48 * t42 * t29 + 0.1e1);
	t25 = (-t28 * t50 - (-t34 * t42 * t43 * t51 + (t37 - 0.1e1) * t46 * t33) * t48 * t54) * t26;
	t1 = [t32, t32, t31, 0, 0; t25, t25, (t47 * t28 - (-t40 * t53 + t34 * t46 + (-t34 * t50 + t53) * t31) * t54) * t41 * t26, 0, 0; t27, t27, -0.1e1 / t41 * t44 * t35 * t50, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 19:56:48
	% EndTime: 2019-12-29 19:56:48
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end