% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPR7
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRPR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(5));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(9));
	t1 = sin(pkin(9));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(5)), 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (35->19), mult. (92->37), div. (0->0), fcn. (112->8), ass. (0->19)
	t6 = sin(pkin(5));
	t9 = sin(qJ(3));
	t20 = t6 * t9;
	t19 = pkin(7) + r_i_i_C(3);
	t11 = cos(qJ(3));
	t18 = t11 * t6;
	t10 = sin(qJ(2));
	t5 = sin(pkin(9));
	t17 = t5 * t10;
	t12 = cos(qJ(2));
	t16 = t5 * t12;
	t7 = cos(pkin(9));
	t15 = t7 * t10;
	t14 = t7 * t12;
	t13 = t11 * r_i_i_C(1) - t9 * r_i_i_C(2) + pkin(2);
	t8 = cos(pkin(5));
	t4 = -t8 * t17 + t14;
	t2 = t8 * t15 + t16;
	t1 = [0, t19 * t4 + t13 * (-t8 * t16 - t15), (t5 * t18 - t4 * t9) * r_i_i_C(1) + (-t4 * t11 - t5 * t20) * r_i_i_C(2), 0, 0; 0, t19 * t2 + t13 * (t8 * t14 - t17), (-t7 * t18 - t2 * t9) * r_i_i_C(1) + (-t2 * t11 + t7 * t20) * r_i_i_C(2), 0, 0; 1, (t19 * t10 + t13 * t12) * t6, (-t10 * t20 + t11 * t8) * r_i_i_C(1) + (-t10 * t18 - t8 * t9) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (84->23), mult. (222->43), div. (0->0), fcn. (279->10), ass. (0->24)
	t15 = sin(pkin(5));
	t19 = sin(qJ(3));
	t30 = t15 * t19;
	t21 = cos(qJ(3));
	t29 = t15 * t21;
	t18 = cos(pkin(5));
	t20 = sin(qJ(2));
	t28 = t18 * t20;
	t22 = cos(qJ(2));
	t27 = t18 * t22;
	t26 = r_i_i_C(3) + qJ(4);
	t13 = sin(pkin(10));
	t16 = cos(pkin(10));
	t25 = r_i_i_C(1) * t16 - r_i_i_C(2) * t13 + pkin(3);
	t24 = t13 * r_i_i_C(1) + t16 * r_i_i_C(2) + pkin(7);
	t23 = t26 * t19 + t25 * t21 + pkin(2);
	t17 = cos(pkin(9));
	t14 = sin(pkin(9));
	t9 = -t18 * t21 + t20 * t30;
	t8 = -t14 * t28 + t17 * t22;
	t6 = t14 * t22 + t17 * t28;
	t3 = -t14 * t29 + t19 * t8;
	t1 = t17 * t29 + t19 * t6;
	t2 = [0, t24 * t8 + t23 * (-t14 * t27 - t17 * t20), t26 * (t14 * t30 + t21 * t8) - t25 * t3, t3, 0; 0, t24 * t6 + t23 * (-t14 * t20 + t17 * t27), t26 * (-t17 * t30 + t21 * t6) - t25 * t1, t1, 0; 1, (t24 * t20 + t23 * t22) * t15, t26 * (t18 * t19 + t20 * t29) - t25 * t9, t9, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (190->49), mult. (510->92), div. (0->0), fcn. (658->12), ass. (0->39)
	t27 = sin(qJ(3));
	t30 = cos(qJ(3));
	t26 = sin(qJ(5));
	t29 = cos(qJ(5));
	t34 = t26 * r_i_i_C(1) + t29 * r_i_i_C(2) + qJ(4);
	t49 = t30 * pkin(3) + t34 * t27 + pkin(2);
	t48 = r_i_i_C(3) + pkin(8);
	t23 = sin(pkin(10));
	t47 = t23 * t30;
	t24 = sin(pkin(5));
	t28 = sin(qJ(2));
	t46 = t24 * t28;
	t25 = cos(pkin(10));
	t45 = t25 * t30;
	t31 = cos(qJ(2));
	t44 = t30 * t31;
	t43 = cos(pkin(5));
	t42 = cos(pkin(9));
	t41 = sin(pkin(9));
	t39 = t24 * t42;
	t38 = t24 * t41;
	t37 = t43 * t42;
	t36 = t43 * t41;
	t35 = t29 * r_i_i_C(1) - t26 * r_i_i_C(2) + pkin(4);
	t32 = -t48 * t23 - t35 * t25 - pkin(3);
	t19 = t43 * t27 + t30 * t46;
	t18 = t27 * t46 - t43 * t30;
	t17 = -t28 * t36 + t42 * t31;
	t16 = t42 * t28 + t31 * t36;
	t15 = t28 * t37 + t41 * t31;
	t14 = t41 * t28 - t31 * t37;
	t11 = t17 * t30 + t27 * t38;
	t10 = t17 * t27 - t30 * t38;
	t9 = t15 * t30 - t27 * t39;
	t8 = t15 * t27 + t30 * t39;
	t7 = -t24 * t31 * t23 + t19 * t25;
	t2 = t11 * t25 + t16 * t23;
	t1 = t14 * t23 + t9 * t25;
	t3 = [0, t17 * pkin(7) + t48 * (-t16 * t47 - t17 * t25) + t35 * (-t16 * t45 + t17 * t23) - t49 * t16, t32 * t10 + t34 * t11, t10, (t10 * t29 - t2 * t26) * r_i_i_C(1) + (-t10 * t26 - t2 * t29) * r_i_i_C(2); 0, t15 * pkin(7) + t48 * (-t14 * t47 - t15 * t25) + t35 * (-t14 * t45 + t15 * t23) - t49 * t14, t32 * t8 + t34 * t9, t8, (-t1 * t26 + t8 * t29) * r_i_i_C(1) + (-t1 * t29 - t8 * t26) * r_i_i_C(2); 1, (t48 * (t23 * t44 - t25 * t28) + t35 * (t23 * t28 + t25 * t44) + t28 * pkin(7) + t49 * t31) * t24, t32 * t18 + t34 * t19, t18, (t18 * t29 - t7 * t26) * r_i_i_C(1) + (-t18 * t26 - t7 * t29) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,5);
end