% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:50
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.15s
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
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (175->29), mult. (265->51), div. (0->0), fcn. (329->10), ass. (0->30)
	t46 = pkin(4) - r_i_i_C(2);
	t45 = r_i_i_C(3) + qJ(5);
	t44 = r_i_i_C(1) + pkin(9) + pkin(8);
	t25 = sin(pkin(11));
	t26 = sin(pkin(6));
	t43 = t25 * t26;
	t27 = cos(pkin(11));
	t42 = t26 * t27;
	t30 = sin(qJ(2));
	t41 = t26 * t30;
	t31 = cos(qJ(3));
	t40 = t26 * t31;
	t28 = cos(pkin(6));
	t39 = t28 * t30;
	t32 = cos(qJ(2));
	t38 = t28 * t32;
	t17 = t25 * t32 + t27 * t39;
	t24 = qJ(3) + qJ(4);
	t22 = sin(t24);
	t23 = cos(t24);
	t7 = t17 * t22 + t23 * t42;
	t37 = t45 * (t17 * t23 - t22 * t42) - t46 * t7;
	t19 = -t25 * t39 + t27 * t32;
	t9 = t19 * t22 - t23 * t43;
	t36 = -t46 * t9 + t45 * (t19 * t23 + t22 * t43);
	t14 = t22 * t41 - t28 * t23;
	t35 = t45 * (t28 * t22 + t23 * t41) - t46 * t14;
	t34 = t31 * pkin(3) + t45 * t22 + t46 * t23 + pkin(2);
	t29 = sin(qJ(3));
	t1 = [0, t44 * t19 + t34 * (-t25 * t38 - t27 * t30), (-t19 * t29 + t25 * t40) * pkin(3) + t36, t36, t9, 0; 0, t44 * t17 + t34 * (-t25 * t30 + t27 * t38), (-t17 * t29 - t27 * t40) * pkin(3) + t37, t37, t7, 0; 1, (t44 * t30 + t34 * t32) * t26, (t28 * t31 - t29 * t41) * pkin(3) + t35, t35, t14, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:50:24
	% EndTime: 2019-10-09 22:50:24
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (292->42), mult. (463->73), div. (0->0), fcn. (581->12), ass. (0->34)
	t59 = pkin(4) + pkin(10) + r_i_i_C(3);
	t33 = sin(qJ(6));
	t36 = cos(qJ(6));
	t58 = t33 * r_i_i_C(1) + t36 * r_i_i_C(2) + qJ(5);
	t30 = qJ(3) + qJ(4);
	t28 = sin(t30);
	t29 = cos(t30);
	t37 = cos(qJ(3));
	t57 = t37 * pkin(3) + t58 * t28 + t59 * t29 + pkin(2);
	t31 = sin(pkin(11));
	t32 = sin(pkin(6));
	t54 = t31 * t32;
	t35 = sin(qJ(2));
	t53 = t32 * t35;
	t38 = cos(qJ(2));
	t52 = t32 * t38;
	t51 = cos(pkin(6));
	t50 = cos(pkin(11));
	t48 = t31 * t51;
	t47 = t32 * t50;
	t46 = t51 * t50;
	t44 = t36 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(5) + pkin(8) + pkin(9);
	t20 = t31 * t38 + t35 * t46;
	t9 = t20 * t28 + t29 * t47;
	t43 = -t59 * t9 + t58 * (t20 * t29 - t28 * t47);
	t22 = -t35 * t48 + t50 * t38;
	t11 = t22 * t28 - t29 * t54;
	t42 = t58 * (t22 * t29 + t28 * t54) - t59 * t11;
	t17 = t28 * t53 - t51 * t29;
	t41 = t58 * (t51 * t28 + t29 * t53) - t59 * t17;
	t34 = sin(qJ(3));
	t21 = t50 * t35 + t38 * t48;
	t19 = t31 * t35 - t38 * t46;
	t1 = [0, -t21 * t57 + t44 * t22, (-t22 * t34 + t37 * t54) * pkin(3) + t42, t42, t11, (t11 * t36 - t21 * t33) * r_i_i_C(1) + (-t11 * t33 - t21 * t36) * r_i_i_C(2); 0, -t19 * t57 + t44 * t20, (-t20 * t34 - t37 * t47) * pkin(3) + t43, t43, t9, (-t19 * t33 + t9 * t36) * r_i_i_C(1) + (-t19 * t36 - t9 * t33) * r_i_i_C(2); 1, (t44 * t35 + t57 * t38) * t32, (-t34 * t53 + t51 * t37) * pkin(3) + t41, t41, t17, (t17 * t36 + t33 * t52) * r_i_i_C(1) + (-t17 * t33 + t36 * t52) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end