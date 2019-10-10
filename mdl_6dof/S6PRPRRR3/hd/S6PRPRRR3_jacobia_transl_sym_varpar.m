% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:57
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (25->11), mult. (63->20), div. (0->0), fcn. (78->8), ass. (0->13)
	t11 = cos(pkin(6));
	t12 = sin(qJ(2));
	t17 = t11 * t12;
	t13 = cos(qJ(2));
	t16 = t11 * t13;
	t15 = r_i_i_C(3) + qJ(3);
	t14 = r_i_i_C(1) * cos(pkin(12)) - r_i_i_C(2) * sin(pkin(12)) + pkin(2);
	t10 = cos(pkin(11));
	t8 = sin(pkin(6));
	t7 = sin(pkin(11));
	t3 = t10 * t12 + t7 * t16;
	t1 = -t10 * t16 + t7 * t12;
	t2 = [0, t15 * (t10 * t13 - t7 * t17) - t14 * t3, t3, 0, 0, 0; 0, t15 * (t10 * t17 + t7 * t13) - t14 * t1, t1, 0, 0, 0; 1, (t15 * t12 + t14 * t13) * t8, -t8 * t13, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (62->22), mult. (102->42), div. (0->0), fcn. (127->9), ass. (0->21)
	t23 = r_i_i_C(3) + pkin(8) + qJ(3);
	t10 = sin(pkin(11));
	t11 = sin(pkin(6));
	t22 = t10 * t11;
	t12 = cos(pkin(11));
	t21 = t11 * t12;
	t15 = sin(qJ(2));
	t20 = t11 * t15;
	t13 = cos(pkin(6));
	t19 = t13 * t15;
	t16 = cos(qJ(2));
	t18 = t13 * t16;
	t9 = pkin(12) + qJ(4);
	t7 = sin(t9);
	t8 = cos(t9);
	t17 = r_i_i_C(1) * t8 - r_i_i_C(2) * t7 + cos(pkin(12)) * pkin(3) + pkin(2);
	t4 = -t10 * t19 + t12 * t16;
	t3 = t10 * t18 + t12 * t15;
	t2 = t10 * t16 + t12 * t19;
	t1 = t10 * t15 - t12 * t18;
	t5 = [0, -t17 * t3 + t23 * t4, t3, (t8 * t22 - t4 * t7) * r_i_i_C(1) + (-t7 * t22 - t4 * t8) * r_i_i_C(2), 0, 0; 0, -t17 * t1 + t23 * t2, t1, (-t2 * t7 - t8 * t21) * r_i_i_C(1) + (-t2 * t8 + t7 * t21) * r_i_i_C(2), 0, 0; 1, (t23 * t15 + t17 * t16) * t11, -t11 * t16, (t13 * t8 - t7 * t20) * r_i_i_C(1) + (-t13 * t7 - t8 * t20) * r_i_i_C(2), 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (140->30), mult. (162->52), div. (0->0), fcn. (199->11), ass. (0->27)
	t19 = pkin(12) + qJ(4);
	t17 = qJ(5) + t19;
	t13 = sin(t17);
	t14 = cos(t17);
	t21 = sin(pkin(6));
	t22 = cos(pkin(11));
	t30 = t21 * t22;
	t20 = sin(pkin(11));
	t25 = cos(qJ(2));
	t23 = cos(pkin(6));
	t24 = sin(qJ(2));
	t28 = t23 * t24;
	t8 = t20 * t25 + t22 * t28;
	t35 = (-t8 * t13 - t14 * t30) * r_i_i_C(1) + (t13 * t30 - t8 * t14) * r_i_i_C(2);
	t10 = -t20 * t28 + t22 * t25;
	t31 = t20 * t21;
	t34 = (-t10 * t13 + t14 * t31) * r_i_i_C(1) + (-t10 * t14 - t13 * t31) * r_i_i_C(2);
	t29 = t21 * t24;
	t33 = (-t13 * t29 + t23 * t14) * r_i_i_C(1) + (-t23 * t13 - t14 * t29) * r_i_i_C(2);
	t32 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3);
	t27 = t23 * t25;
	t16 = cos(t19);
	t26 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + pkin(4) * t16 + cos(pkin(12)) * pkin(3) + pkin(2);
	t15 = sin(t19);
	t9 = t20 * t27 + t22 * t24;
	t7 = t20 * t24 - t22 * t27;
	t1 = [0, t32 * t10 - t26 * t9, t9, (-t10 * t15 + t16 * t31) * pkin(4) + t34, t34, 0; 0, -t26 * t7 + t32 * t8, t7, (-t15 * t8 - t16 * t30) * pkin(4) + t35, t35, 0; 1, (t32 * t24 + t26 * t25) * t21, -t21 * t25, (-t15 * t29 + t16 * t23) * pkin(4) + t33, t33, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (345->45), mult. (399->72), div. (0->0), fcn. (500->13), ass. (0->37)
	t57 = pkin(10) + r_i_i_C(3);
	t33 = sin(qJ(6));
	t35 = cos(qJ(6));
	t56 = t35 * r_i_i_C(1) - t33 * r_i_i_C(2) + pkin(5);
	t29 = pkin(12) + qJ(4);
	t27 = qJ(5) + t29;
	t23 = sin(t27);
	t24 = cos(t27);
	t26 = cos(t29);
	t55 = t57 * t23 + t56 * t24 + pkin(4) * t26 + cos(pkin(12)) * pkin(3) + pkin(2);
	t30 = sin(pkin(11));
	t31 = sin(pkin(6));
	t51 = t30 * t31;
	t34 = sin(qJ(2));
	t50 = t30 * t34;
	t36 = cos(qJ(2));
	t49 = t30 * t36;
	t48 = t31 * t34;
	t47 = t31 * t36;
	t46 = cos(pkin(11));
	t45 = t31 * t46;
	t44 = t46 * t34;
	t43 = t46 * t36;
	t41 = t33 * r_i_i_C(1) + t35 * r_i_i_C(2) + pkin(8) + pkin(9) + qJ(3);
	t32 = cos(pkin(6));
	t17 = t32 * t44 + t49;
	t8 = t17 * t24 - t23 * t45;
	t40 = t57 * t8 + t56 * (-t17 * t23 - t24 * t45);
	t19 = -t32 * t50 + t43;
	t10 = t19 * t24 + t23 * t51;
	t39 = t57 * t10 + t56 * (-t19 * t23 + t24 * t51);
	t15 = t32 * t23 + t24 * t48;
	t38 = t57 * t15 + t56 * (-t23 * t48 + t32 * t24);
	t25 = sin(t29);
	t18 = t32 * t49 + t44;
	t16 = -t32 * t43 + t50;
	t1 = [0, -t18 * t55 + t41 * t19, t18, (-t19 * t25 + t26 * t51) * pkin(4) + t39, t39, (-t10 * t33 + t18 * t35) * r_i_i_C(1) + (-t10 * t35 - t18 * t33) * r_i_i_C(2); 0, -t16 * t55 + t41 * t17, t16, (-t17 * t25 - t26 * t45) * pkin(4) + t40, t40, (t16 * t35 - t8 * t33) * r_i_i_C(1) + (-t16 * t33 - t8 * t35) * r_i_i_C(2); 1, (t41 * t34 + t55 * t36) * t31, -t47, (-t25 * t48 + t26 * t32) * pkin(4) + t38, t38, (-t15 * t33 - t35 * t47) * r_i_i_C(1) + (-t15 * t35 + t33 * t47) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end