% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:46
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:42
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:42
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:42
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(11));
	t1 = sin(pkin(11));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:42
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (35->19), mult. (92->37), div. (0->0), fcn. (112->8), ass. (0->19)
	t6 = sin(pkin(6));
	t9 = sin(qJ(3));
	t20 = t6 * t9;
	t19 = pkin(8) + r_i_i_C(3);
	t11 = cos(qJ(3));
	t18 = t11 * t6;
	t10 = sin(qJ(2));
	t5 = sin(pkin(11));
	t17 = t5 * t10;
	t12 = cos(qJ(2));
	t16 = t5 * t12;
	t7 = cos(pkin(11));
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
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:42
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (92->27), mult. (152->51), div. (0->0), fcn. (184->10), ass. (0->25)
	t14 = qJ(3) + qJ(4);
	t12 = sin(t14);
	t13 = cos(t14);
	t16 = sin(pkin(6));
	t17 = cos(pkin(11));
	t29 = t16 * t17;
	t15 = sin(pkin(11));
	t22 = cos(qJ(2));
	t18 = cos(pkin(6));
	t20 = sin(qJ(2));
	t26 = t18 * t20;
	t8 = t15 * t22 + t17 * t26;
	t34 = (-t12 * t8 - t13 * t29) * r_i_i_C(1) + (t12 * t29 - t13 * t8) * r_i_i_C(2);
	t10 = -t15 * t26 + t17 * t22;
	t30 = t15 * t16;
	t33 = (-t10 * t12 + t13 * t30) * r_i_i_C(1) + (-t10 * t13 - t12 * t30) * r_i_i_C(2);
	t28 = t16 * t20;
	t32 = (-t12 * t28 + t13 * t18) * r_i_i_C(1) + (-t12 * t18 - t13 * t28) * r_i_i_C(2);
	t31 = r_i_i_C(3) + pkin(9) + pkin(8);
	t21 = cos(qJ(3));
	t27 = t16 * t21;
	t25 = t18 * t22;
	t24 = pkin(3) * t21 + r_i_i_C(1) * t13 - r_i_i_C(2) * t12 + pkin(2);
	t19 = sin(qJ(3));
	t1 = [0, t31 * t10 + t24 * (-t15 * t25 - t17 * t20), (-t10 * t19 + t15 * t27) * pkin(3) + t33, t33, 0, 0; 0, t31 * t8 + t24 * (-t15 * t20 + t17 * t25), (-t17 * t27 - t19 * t8) * pkin(3) + t34, t34, 0, 0; 1, (t31 * t20 + t24 * t22) * t16, (t18 * t21 - t19 * t28) * pkin(3) + t32, t32, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (160->37), mult. (189->60), div. (0->0), fcn. (228->12), ass. (0->29)
	t23 = qJ(3) + qJ(4);
	t18 = pkin(12) + t23;
	t15 = sin(t18);
	t16 = cos(t18);
	t25 = sin(pkin(6));
	t26 = cos(pkin(11));
	t34 = t25 * t26;
	t24 = sin(pkin(11));
	t29 = cos(qJ(2));
	t27 = cos(pkin(6));
	t28 = sin(qJ(2));
	t32 = t27 * t28;
	t8 = t24 * t29 + t26 * t32;
	t39 = (-t15 * t8 - t16 * t34) * r_i_i_C(1) + (t15 * t34 - t16 * t8) * r_i_i_C(2);
	t10 = -t24 * t32 + t26 * t29;
	t35 = t24 * t25;
	t38 = (-t10 * t15 + t16 * t35) * r_i_i_C(1) + (-t10 * t16 - t15 * t35) * r_i_i_C(2);
	t33 = t25 * t28;
	t37 = (-t15 * t33 + t16 * t27) * r_i_i_C(1) + (-t15 * t27 - t16 * t33) * r_i_i_C(2);
	t36 = r_i_i_C(3) + qJ(5) + pkin(9) + pkin(8);
	t31 = t27 * t29;
	t20 = cos(t23);
	t13 = pkin(4) * t20 + cos(qJ(3)) * pkin(3);
	t30 = r_i_i_C(1) * t16 - r_i_i_C(2) * t15 + pkin(2) + t13;
	t19 = sin(t23);
	t12 = -pkin(4) * t19 - sin(qJ(3)) * pkin(3);
	t9 = t24 * t31 + t26 * t28;
	t7 = t24 * t28 - t26 * t31;
	t1 = [0, t36 * t10 - t30 * t9, t10 * t12 + t13 * t35 + t38, (-t10 * t19 + t20 * t35) * pkin(4) + t38, t9, 0; 0, -t30 * t7 + t36 * t8, t12 * t8 - t13 * t34 + t39, (-t19 * t8 - t20 * t34) * pkin(4) + t39, t7, 0; 1, (t36 * t28 + t30 * t29) * t25, t12 * t33 + t13 * t27 + t37, (-t19 * t33 + t20 * t27) * pkin(4) + t37, -t25 * t29, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:46:42
	% EndTime: 2019-10-09 22:46:43
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (365->52), mult. (426->82), div. (0->0), fcn. (529->14), ass. (0->37)
	t59 = pkin(10) + r_i_i_C(3);
	t38 = sin(qJ(6));
	t40 = cos(qJ(6));
	t58 = t40 * r_i_i_C(1) - t38 * r_i_i_C(2) + pkin(5);
	t33 = qJ(3) + qJ(4);
	t30 = cos(t33);
	t23 = pkin(4) * t30 + cos(qJ(3)) * pkin(3);
	t28 = pkin(12) + t33;
	t25 = sin(t28);
	t26 = cos(t28);
	t57 = t59 * t25 + t58 * t26 + pkin(2) + t23;
	t34 = sin(pkin(11));
	t35 = sin(pkin(6));
	t53 = t34 * t35;
	t36 = cos(pkin(11));
	t52 = t35 * t36;
	t39 = sin(qJ(2));
	t51 = t35 * t39;
	t41 = cos(qJ(2));
	t50 = t35 * t41;
	t37 = cos(pkin(6));
	t49 = t37 * t39;
	t48 = t37 * t41;
	t46 = t38 * r_i_i_C(1) + t40 * r_i_i_C(2) + pkin(8) + pkin(9) + qJ(5);
	t17 = t34 * t41 + t36 * t49;
	t8 = t17 * t26 - t25 * t52;
	t45 = t59 * t8 + t58 * (-t17 * t25 - t26 * t52);
	t19 = -t34 * t49 + t36 * t41;
	t10 = t19 * t26 + t25 * t53;
	t44 = t59 * t10 + t58 * (-t19 * t25 + t26 * t53);
	t15 = t37 * t25 + t26 * t51;
	t43 = t59 * t15 + t58 * (-t25 * t51 + t37 * t26);
	t29 = sin(t33);
	t22 = -pkin(4) * t29 - sin(qJ(3)) * pkin(3);
	t18 = t34 * t48 + t36 * t39;
	t16 = t34 * t39 - t36 * t48;
	t1 = [0, -t18 * t57 + t46 * t19, t19 * t22 + t23 * t53 + t44, (-t19 * t29 + t30 * t53) * pkin(4) + t44, t18, (-t10 * t38 + t18 * t40) * r_i_i_C(1) + (-t10 * t40 - t18 * t38) * r_i_i_C(2); 0, -t16 * t57 + t46 * t17, t17 * t22 - t23 * t52 + t45, (-t17 * t29 - t30 * t52) * pkin(4) + t45, t16, (t16 * t40 - t8 * t38) * r_i_i_C(1) + (-t16 * t38 - t8 * t40) * r_i_i_C(2); 1, (t46 * t39 + t57 * t41) * t35, t22 * t51 + t37 * t23 + t43, (-t29 * t51 + t30 * t37) * pkin(4) + t43, -t50, (-t15 * t38 - t40 * t50) * r_i_i_C(1) + (-t15 * t40 + t38 * t50) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end