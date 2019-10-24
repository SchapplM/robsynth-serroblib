% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPRP6
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
%   Wie in S5PRPRP6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:25
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRPRP6_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP6_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_jacobia_rot_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:11
	% EndTime: 2019-10-24 10:25:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:11
	% EndTime: 2019-10-24 10:25:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:11
	% EndTime: 2019-10-24 10:25:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:11
	% EndTime: 2019-10-24 10:25:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (45->12), mult. (82->22), div. (24->8), fcn. (122->7), ass. (0->18)
	t22 = cos(qJ(2));
	t19 = sin(pkin(7));
	t21 = sin(qJ(2));
	t25 = t19 * t21;
	t14 = atan2(-t25, -t22);
	t12 = sin(t14);
	t27 = t12 * t22;
	t17 = t21 ^ 2;
	t24 = t22 ^ 2;
	t26 = t17 / t24;
	t23 = t19 ^ 2;
	t20 = cos(pkin(7));
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
	% StartTime: 2019-10-24 10:25:11
	% EndTime: 2019-10-24 10:25:11
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (51->13), mult. (135->34), div. (32->8), fcn. (192->9), ass. (0->23)
	t29 = sin(qJ(2));
	t26 = sin(pkin(7));
	t31 = cos(qJ(2));
	t34 = t26 * t31;
	t22 = atan2(-t34, t29);
	t20 = sin(t22);
	t36 = t20 * t29;
	t25 = t31 ^ 2;
	t35 = 0.1e1 / t29 ^ 2 * t25;
	t27 = cos(pkin(7));
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
	% StartTime: 2019-10-24 10:25:11
	% EndTime: 2019-10-24 10:25:11
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (112->20), mult. (359->57), div. (75->11), fcn. (540->9), ass. (0->34)
	t46 = cos(qJ(2));
	t55 = t46 ^ 2;
	t41 = sin(pkin(7));
	t42 = cos(pkin(7));
	t43 = sin(qJ(4));
	t44 = sin(qJ(2));
	t45 = cos(qJ(4));
	t49 = t44 * t45;
	t33 = t41 * t49 + t42 * t43;
	t48 = t46 * t45;
	t29 = atan2(t33, t48);
	t25 = sin(t29);
	t26 = cos(t29);
	t24 = t25 * t33 + t26 * t48;
	t23 = 0.1e1 / t24 ^ 2;
	t31 = t41 * t43 - t42 * t49;
	t54 = t23 * t31;
	t52 = t26 * t33;
	t51 = t42 * t46;
	t50 = t43 * t44;
	t32 = t41 * t45 + t42 * t50;
	t30 = 0.1e1 / t32 ^ 2;
	t47 = t42 ^ 2 * t55 * t30;
	t40 = 0.1e1 / t55;
	t37 = 0.1e1 / t45 ^ 2;
	t36 = 0.1e1 / t45;
	t34 = -t41 * t50 + t42 * t45;
	t28 = 0.1e1 / (t33 ^ 2 * t40 * t37 + 0.1e1);
	t27 = 0.1e1 / (0.1e1 + t47);
	t22 = 0.1e1 / t24;
	t21 = (t33 * t36 * t40 * t44 + t41) * t28;
	t20 = 0.1e1 / (t31 ^ 2 * t23 + 0.1e1);
	t19 = (t33 * t37 * t43 + t34 * t36) / t46 * t28;
	t1 = [0, t21, 0, t19, 0; 0, (-t21 * t52 * t54 + (-t22 * t51 - (-t26 * t44 + (-t21 + t41) * t46 * t25) * t54) * t45) * t20, 0, (t32 * t22 - (-t26 * t46 * t43 + t25 * t34 + (-t25 * t48 + t52) * t19) * t54) * t20, 0; 0, (t42 * t44 / t32 + t43 * t47) * t27, 0, -t31 * t30 * t27 * t51, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end