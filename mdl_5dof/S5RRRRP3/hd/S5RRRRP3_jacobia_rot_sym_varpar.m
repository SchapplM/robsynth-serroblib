% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRRP3
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
%   Wie in S5RRRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRRRP3_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP3_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->0), mult. (12->0), div. (10->0), fcn. (12->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (45->0), mult. (18->0), div. (15->0), fcn. (18->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 1, 1, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (430->17), mult. (271->39), div. (95->9), fcn. (438->7), ass. (0->30)
	t47 = qJ(1) + qJ(2) + qJ(3);
	t46 = cos(t47);
	t53 = t46 ^ 2;
	t52 = cos(qJ(4));
	t45 = sin(t47);
	t51 = sin(qJ(4));
	t55 = t45 * t51;
	t41 = atan2(-t55, -t52);
	t39 = sin(t41);
	t40 = cos(t41);
	t35 = -t39 * t55 - t40 * t52;
	t34 = 0.1e1 / t35 ^ 2;
	t59 = t34 * t51;
	t58 = t39 * t52;
	t43 = t45 ^ 2;
	t57 = t43 / t53;
	t48 = t51 ^ 2;
	t50 = 0.1e1 / t52 ^ 2;
	t54 = t48 * t50;
	t42 = 0.1e1 / (t43 * t54 + 0.1e1);
	t56 = t45 * t42;
	t49 = 0.1e1 / t52;
	t38 = 0.1e1 / (t50 * t57 + 0.1e1);
	t37 = t46 * t51 * t49 * t42;
	t36 = (0.1e1 + t54) * t56;
	t33 = 0.1e1 / t35;
	t32 = (-0.1e1 - t57) * t49 * t38;
	t31 = 0.1e1 / (t53 * t48 * t34 + 0.1e1);
	t30 = (-t33 * t55 - (-t40 * t48 * t49 * t56 + (t42 - 0.1e1) * t51 * t39) * t53 * t59) * t31;
	t1 = [t37, t37, t37, t36, 0; t30, t30, t30, (t52 * t33 - (-t45 * t58 + t40 * t51 + (-t40 * t55 + t58) * t36) * t59) * t46 * t31, 0; t32, t32, t32, -0.1e1 / t46 * t50 * t38 * t55, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end