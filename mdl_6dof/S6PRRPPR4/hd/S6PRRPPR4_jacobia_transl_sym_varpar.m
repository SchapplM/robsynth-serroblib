% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPPR4
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:43
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:43
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:43
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(10));
	t1 = sin(pkin(10));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:43
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (35->19), mult. (92->37), div. (0->0), fcn. (112->8), ass. (0->19)
	t6 = sin(pkin(6));
	t9 = sin(qJ(3));
	t20 = t6 * t9;
	t19 = pkin(8) + r_i_i_C(3);
	t11 = cos(qJ(3));
	t18 = t11 * t6;
	t10 = sin(qJ(2));
	t5 = sin(pkin(10));
	t17 = t5 * t10;
	t12 = cos(qJ(2));
	t16 = t5 * t12;
	t7 = cos(pkin(10));
	t15 = t7 * t10;
	t14 = t7 * t12;
	t13 = t11 * r_i_i_C(1) - t9 * r_i_i_C(2) + pkin(2);
	t8 = cos(pkin(6));
	t4 = -t8 * t17 + t14;
	t2 = t8 * t15 + t16;
	t1 = [0, t19 * t4 + t13 * (-t8 * t16 - t15), (t5 * t18 - t4 * t9) * r_i_i_C(1) + (-t11 * t4 - t5 * t20) * r_i_i_C(2), 0, 0, 0; 0, t19 * t2 + t13 * (t8 * t14 - t17), (-t7 * t18 - t2 * t9) * r_i_i_C(1) + (-t11 * t2 + t7 * t20) * r_i_i_C(2), 0, 0, 0; 1, (t19 * t10 + t13 * t12) * t6, (-t10 * t20 + t11 * t8) * r_i_i_C(1) + (-t10 * t18 - t8 * t9) * r_i_i_C(2), 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:43
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (84->23), mult. (222->43), div. (0->0), fcn. (279->10), ass. (0->24)
	t15 = sin(pkin(6));
	t19 = sin(qJ(3));
	t30 = t15 * t19;
	t21 = cos(qJ(3));
	t29 = t15 * t21;
	t18 = cos(pkin(6));
	t20 = sin(qJ(2));
	t28 = t18 * t20;
	t22 = cos(qJ(2));
	t27 = t18 * t22;
	t26 = r_i_i_C(3) + qJ(4);
	t13 = sin(pkin(11));
	t16 = cos(pkin(11));
	t25 = r_i_i_C(1) * t16 - r_i_i_C(2) * t13 + pkin(3);
	t24 = t13 * r_i_i_C(1) + t16 * r_i_i_C(2) + pkin(8);
	t23 = t26 * t19 + t25 * t21 + pkin(2);
	t17 = cos(pkin(10));
	t14 = sin(pkin(10));
	t9 = -t18 * t21 + t20 * t30;
	t8 = -t14 * t28 + t17 * t22;
	t6 = t14 * t22 + t17 * t28;
	t3 = -t14 * t29 + t19 * t8;
	t1 = t17 * t29 + t19 * t6;
	t2 = [0, t24 * t8 + t23 * (-t14 * t27 - t17 * t20), t26 * (t14 * t30 + t21 * t8) - t25 * t3, t3, 0, 0; 0, t24 * t6 + t23 * (-t14 * t20 + t17 * t27), t26 * (-t17 * t30 + t21 * t6) - t25 * t1, t1, 0, 0; 1, (t24 * t20 + t23 * t22) * t15, t26 * (t18 * t19 + t20 * t29) - t25 * t9, t9, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:43
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (130->38), mult. (344->70), div. (0->0), fcn. (438->10), ass. (0->34)
	t24 = sin(qJ(3));
	t26 = cos(qJ(3));
	t37 = r_i_i_C(2) + qJ(4);
	t44 = pkin(3) * t26 + t37 * t24 + pkin(2);
	t43 = pkin(4) + r_i_i_C(1);
	t20 = sin(pkin(11));
	t42 = t20 * t26;
	t22 = sin(pkin(6));
	t41 = t22 * t24;
	t40 = t22 * t26;
	t23 = cos(pkin(11));
	t39 = t23 * t26;
	t27 = cos(qJ(2));
	t38 = t26 * t27;
	t36 = r_i_i_C(3) + qJ(5);
	t35 = cos(pkin(6));
	t34 = cos(pkin(10));
	t21 = sin(pkin(10));
	t32 = t21 * t35;
	t31 = t22 * t34;
	t30 = t35 * t34;
	t28 = -t36 * t20 - t43 * t23 - pkin(3);
	t25 = sin(qJ(2));
	t16 = t35 * t24 + t25 * t40;
	t15 = t25 * t41 - t35 * t26;
	t14 = -t25 * t32 + t34 * t27;
	t13 = t34 * t25 + t27 * t32;
	t12 = t21 * t27 + t25 * t30;
	t11 = t21 * t25 - t27 * t30;
	t8 = t14 * t26 + t21 * t41;
	t7 = t14 * t24 - t21 * t40;
	t6 = t12 * t26 - t24 * t31;
	t5 = t12 * t24 + t26 * t31;
	t1 = [0, t14 * pkin(8) + t43 * (-t13 * t39 + t14 * t20) + t36 * (-t13 * t42 - t14 * t23) - t44 * t13, t28 * t7 + t37 * t8, t7, -t13 * t23 + t8 * t20, 0; 0, t12 * pkin(8) + t43 * (-t11 * t39 + t12 * t20) + t36 * (-t11 * t42 - t12 * t23) - t44 * t11, t28 * t5 + t37 * t6, t5, -t11 * t23 + t6 * t20, 0; 1, (t36 * (t20 * t38 - t23 * t25) + t43 * (t20 * t25 + t23 * t38) + pkin(8) * t25 + t44 * t27) * t22, t28 * t15 + t37 * t16, t15, t22 * t27 * t23 + t16 * t20, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:12:43
	% EndTime: 2019-10-09 22:12:44
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (246->53), mult. (658->98), div. (0->0), fcn. (853->12), ass. (0->43)
	t31 = sin(qJ(3));
	t47 = -r_i_i_C(3) - pkin(9) + qJ(4);
	t53 = cos(qJ(3));
	t54 = pkin(3) * t53 + t31 * t47 + pkin(2);
	t28 = sin(pkin(6));
	t32 = sin(qJ(2));
	t52 = t28 * t32;
	t34 = cos(qJ(2));
	t51 = t28 * t34;
	t50 = cos(pkin(6));
	t49 = cos(pkin(10));
	t48 = sin(pkin(10));
	t27 = sin(pkin(11));
	t46 = t27 * t53;
	t29 = cos(pkin(11));
	t45 = t29 * t53;
	t44 = t34 * t53;
	t43 = t28 * t49;
	t42 = t28 * t48;
	t40 = t50 * t49;
	t39 = t50 * t48;
	t30 = sin(qJ(6));
	t33 = cos(qJ(6));
	t38 = r_i_i_C(1) * t30 + r_i_i_C(2) * t33 + qJ(5);
	t37 = r_i_i_C(1) * t33 - r_i_i_C(2) * t30 + pkin(4) + pkin(5);
	t35 = -t27 * t38 - t29 * t37 - pkin(3);
	t22 = t31 * t50 + t52 * t53;
	t21 = t31 * t52 - t50 * t53;
	t20 = -t32 * t39 + t34 * t49;
	t19 = t32 * t49 + t34 * t39;
	t18 = t32 * t40 + t34 * t48;
	t17 = t32 * t48 - t34 * t40;
	t14 = t20 * t53 + t31 * t42;
	t13 = t20 * t31 - t42 * t53;
	t12 = t18 * t53 - t31 * t43;
	t11 = t18 * t31 + t43 * t53;
	t10 = t22 * t29 - t27 * t51;
	t9 = t22 * t27 + t29 * t51;
	t4 = t14 * t29 + t19 * t27;
	t3 = t14 * t27 - t19 * t29;
	t2 = t12 * t29 + t17 * t27;
	t1 = t12 * t27 - t17 * t29;
	t5 = [0, t20 * pkin(8) + t38 * (-t19 * t46 - t20 * t29) + t37 * (-t19 * t45 + t20 * t27) - t54 * t19, t13 * t35 + t14 * t47, t13, t3, (t3 * t33 - t30 * t4) * r_i_i_C(1) + (-t3 * t30 - t33 * t4) * r_i_i_C(2); 0, t18 * pkin(8) + t38 * (-t17 * t46 - t18 * t29) + t37 * (-t17 * t45 + t18 * t27) - t54 * t17, t11 * t35 + t12 * t47, t11, t1, (t1 * t33 - t2 * t30) * r_i_i_C(1) + (-t1 * t30 - t2 * t33) * r_i_i_C(2); 1, (t38 * (t27 * t44 - t29 * t32) + t37 * (t27 * t32 + t29 * t44) + t32 * pkin(8) + t54 * t34) * t28, t21 * t35 + t22 * t47, t21, t9, (-t10 * t30 + t33 * t9) * r_i_i_C(1) + (-t10 * t33 - t30 * t9) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end