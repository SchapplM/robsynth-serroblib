% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PPPRR1
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
%   Wie in S5PPPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PPPRR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:23
	% EndTime: 2019-12-05 14:58:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:23
	% EndTime: 2019-12-05 14:58:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (24->0), mult. (27->0), div. (5->0), fcn. (35->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 1, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 14:58:24
	% EndTime: 2019-12-05 14:58:24
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (222->18), mult. (268->42), div. (47->11), fcn. (395->11), ass. (0->31)
	t42 = pkin(9) + qJ(4);
	t39 = sin(t42);
	t43 = sin(pkin(8));
	t54 = t43 * t39;
	t46 = cos(pkin(7));
	t53 = t43 * t46;
	t44 = sin(pkin(7));
	t45 = cos(pkin(8));
	t52 = t44 * t45;
	t51 = t46 * t39;
	t40 = cos(t42);
	t50 = t46 * t40;
	t36 = t44 * t39 + t45 * t50;
	t47 = sin(qJ(5));
	t48 = cos(qJ(5));
	t27 = t36 * t48 + t47 * t53;
	t25 = 0.1e1 / t27 ^ 2;
	t26 = t36 * t47 - t48 * t53;
	t49 = t26 ^ 2 * t25 + 0.1e1;
	t38 = 0.1e1 / t39 ^ 2;
	t35 = -t44 * t40 + t45 * t51;
	t34 = t40 * t52 - t51;
	t32 = t39 * t52 + t50;
	t31 = atan2(-t32, t54);
	t29 = cos(t31);
	t28 = sin(t31);
	t24 = 0.1e1 / t49;
	t23 = -t28 * t32 + t29 * t54;
	t22 = 0.1e1 / t23 ^ 2;
	t20 = (-t34 / t39 + t40 * t32 * t38) / t43 / (0.1e1 + t32 ^ 2 / t43 ^ 2 * t38);
	t1 = [0, 0, 0, t20, 0; 0, 0, 0, (t36 / t23 - (t29 * t43 * t40 - t28 * t34 + (-t28 * t54 - t29 * t32) * t20) * t35 * t22) / (t35 ^ 2 * t22 + 0.1e1), 0; 0, 0, 0, (-t47 / t27 + t48 * t26 * t25) * t35 * t24, t49 * t24;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end