% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPRR7
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
%   Wie in S5PRPRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRPRR7_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR7_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (45->12), mult. (82->22), div. (24->8), fcn. (122->7), ass. (0->18)
	t22 = cos(qJ(2));
	t19 = sin(pkin(8));
	t21 = sin(qJ(2));
	t25 = t19 * t21;
	t14 = atan2(-t25, -t22);
	t12 = sin(t14);
	t27 = t12 * t22;
	t17 = t21 ^ 2;
	t24 = t22 ^ 2;
	t26 = t17 / t24;
	t23 = t19 ^ 2;
	t20 = cos(pkin(8));
	t16 = t20 ^ 2;
	t13 = cos(t14);
	t11 = (0.1e1 + t26) * t19 / (t23 * t26 + 0.1e1);
	t10 = -t12 * t25 - t13 * t22;
	t9 = 0.1e1 / t10 ^ 2;
	t1 = [0, t11, 0, 0, 0; 0, (t22 / t10 - (-t19 * t27 + t13 * t21 + (-t13 * t25 + t27) * t11) * t21 * t9) / (t16 * t17 * t9 + 0.1e1) * t20, 0, 0, 0; 0, -t20 * t21 / t19 / (0.1e1 + t16 * t24 / t23), 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:10
	% EndTime: 2019-12-05 16:01:11
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (51->13), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->23)
	t29 = sin(qJ(2));
	t26 = sin(pkin(8));
	t31 = cos(qJ(2));
	t34 = t26 * t31;
	t22 = atan2(-t34, t29);
	t20 = sin(t22);
	t36 = t20 * t29;
	t25 = t31 ^ 2;
	t35 = 0.1e1 / t29 ^ 2 * t25;
	t27 = cos(pkin(8));
	t33 = t27 * t29;
	t28 = sin(qJ(4));
	t30 = cos(qJ(4));
	t19 = t26 * t30 + t28 * t33;
	t17 = 0.1e1 / t19 ^ 2;
	t18 = t26 * t28 - t30 * t33;
	t32 = t18 ^ 2 * t17 + 0.1e1;
	t21 = cos(t22);
	t16 = (0.1e1 + t35) * t26 / (t26 ^ 2 * t35 + 0.1e1);
	t15 = 0.1e1 / t32;
	t14 = -t20 * t34 + t21 * t29;
	t13 = 0.1e1 / t14 ^ 2;
	t1 = [0, t16, 0, 0, 0; 0, (-t29 / t14 - (t26 * t36 + t21 * t31 + (-t21 * t34 - t36) * t16) * t31 * t13) * t27 / (t27 ^ 2 * t25 * t13 + 0.1e1), 0, 0, 0; 0, (-t30 / t19 - t28 * t18 * t17) * t31 * t27 * t15, 0, t32 * t15, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:01:11
	% EndTime: 2019-12-05 16:01:11
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (105->14), mult. (162->34), div. (37->8), fcn. (227->9), ass. (0->25)
	t45 = sin(qJ(2));
	t43 = sin(pkin(8));
	t46 = cos(qJ(2));
	t49 = t43 * t46;
	t36 = atan2(-t49, t45);
	t34 = sin(t36);
	t51 = t34 * t45;
	t41 = t46 ^ 2;
	t50 = 0.1e1 / t45 ^ 2 * t41;
	t44 = cos(pkin(8));
	t48 = t44 * t45;
	t42 = qJ(4) + qJ(5);
	t38 = sin(t42);
	t39 = cos(t42);
	t33 = t38 * t48 + t43 * t39;
	t31 = 0.1e1 / t33 ^ 2;
	t32 = t43 * t38 - t39 * t48;
	t47 = t32 ^ 2 * t31 + 0.1e1;
	t35 = cos(t36);
	t30 = (0.1e1 + t50) * t43 / (t43 ^ 2 * t50 + 0.1e1);
	t29 = -t34 * t49 + t35 * t45;
	t28 = 0.1e1 / t29 ^ 2;
	t27 = 0.1e1 / t47;
	t25 = t47 * t27;
	t1 = [0, t30, 0, 0, 0; 0, (-t45 / t29 - (t43 * t51 + t35 * t46 + (-t35 * t49 - t51) * t30) * t46 * t28) * t44 / (t44 ^ 2 * t41 * t28 + 0.1e1), 0, 0, 0; 0, (-t39 / t33 - t38 * t32 * t31) * t46 * t44 * t27, 0, t25, t25;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end