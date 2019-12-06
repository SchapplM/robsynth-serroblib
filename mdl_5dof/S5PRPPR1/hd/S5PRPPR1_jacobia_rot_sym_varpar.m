% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPPR1
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
%   Wie in S5PRPPR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRPPR1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR1_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_jacobia_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:22:46
	% EndTime: 2019-12-05 15:22:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:22:45
	% EndTime: 2019-12-05 15:22:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:22:46
	% EndTime: 2019-12-05 15:22:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:22:45
	% EndTime: 2019-12-05 15:22:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:22:46
	% EndTime: 2019-12-05 15:22:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (91->15), mult. (89->34), div. (20->9), fcn. (140->9), ass. (0->22)
	t31 = cos(pkin(8));
	t27 = pkin(7) + qJ(2);
	t23 = sin(t27);
	t29 = sin(pkin(8));
	t34 = t23 * t29;
	t21 = atan2(-t34, -t31);
	t19 = sin(t21);
	t20 = cos(t21);
	t14 = -t19 * t34 - t20 * t31;
	t24 = cos(t27);
	t36 = 0.1e1 / t14 ^ 2 * t24 ^ 2;
	t25 = t29 ^ 2;
	t22 = 0.1e1 / (0.1e1 + t23 ^ 2 * t25 / t31 ^ 2);
	t35 = t22 / t31;
	t28 = sin(pkin(9));
	t33 = t28 * t31;
	t30 = cos(pkin(9));
	t32 = t30 * t31;
	t18 = t23 * t28 + t24 * t32;
	t17 = -t23 * t30 + t24 * t33;
	t16 = 0.1e1 / t18 ^ 2;
	t1 = [0, t24 * t29 * t35, 0, 0, 0; 0, (-0.1e1 / t14 * t34 - (-t20 * t23 * t25 * t35 + (t22 - 0.1e1) * t29 * t19) * t29 * t36) / (t25 * t36 + 0.1e1), 0, 0, 0; 0, ((-t23 * t33 - t24 * t30) / t18 - (-t23 * t32 + t24 * t28) * t17 * t16) / (t17 ^ 2 * t16 + 0.1e1), 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:22:46
	% EndTime: 2019-12-05 15:22:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (147->16), mult. (116->36), div. (25->9), fcn. (175->9), ass. (0->25)
	t42 = cos(pkin(8));
	t40 = pkin(7) + qJ(2);
	t34 = sin(t40);
	t41 = sin(pkin(8));
	t46 = t34 * t41;
	t31 = atan2(-t46, -t42);
	t29 = sin(t31);
	t30 = cos(t31);
	t25 = -t29 * t46 - t30 * t42;
	t36 = cos(t40);
	t48 = 0.1e1 / t25 ^ 2 * t36 ^ 2;
	t37 = t41 ^ 2;
	t32 = 0.1e1 / (0.1e1 + t34 ^ 2 * t37 / t42 ^ 2);
	t47 = t32 / t42;
	t45 = t34 * t42;
	t44 = t36 * t42;
	t39 = pkin(9) + qJ(5);
	t33 = sin(t39);
	t35 = cos(t39);
	t28 = t34 * t33 + t35 * t44;
	t26 = 0.1e1 / t28 ^ 2;
	t27 = t33 * t44 - t34 * t35;
	t43 = t27 ^ 2 * t26 + 0.1e1;
	t23 = 0.1e1 / t43;
	t1 = [0, t36 * t41 * t47, 0, 0, 0; 0, (-0.1e1 / t25 * t46 - (-t30 * t34 * t37 * t47 + (t32 - 0.1e1) * t41 * t29) * t41 * t48) / (t37 * t48 + 0.1e1), 0, 0, 0; 0, ((-t33 * t45 - t36 * t35) / t28 - (t36 * t33 - t35 * t45) * t27 * t26) * t23, 0, 0, t43 * t23;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end