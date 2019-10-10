% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPPR2
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
% Datum: 2019-10-09 22:09
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
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
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
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
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (70->28), mult. (123->52), div. (0->0), fcn. (150->10), ass. (0->24)
	t26 = r_i_i_C(3) + qJ(4) + pkin(8);
	t10 = sin(pkin(10));
	t11 = sin(pkin(6));
	t25 = t10 * t11;
	t12 = cos(pkin(10));
	t24 = t11 * t12;
	t16 = sin(qJ(2));
	t23 = t11 * t16;
	t17 = cos(qJ(3));
	t22 = t11 * t17;
	t13 = cos(pkin(6));
	t21 = t13 * t16;
	t18 = cos(qJ(2));
	t20 = t13 * t18;
	t9 = qJ(3) + pkin(11);
	t7 = sin(t9);
	t8 = cos(t9);
	t19 = t17 * pkin(3) + t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + pkin(2);
	t15 = sin(qJ(3));
	t4 = -t10 * t21 + t12 * t18;
	t3 = t10 * t20 + t12 * t16;
	t2 = t10 * t18 + t12 * t21;
	t1 = t10 * t16 - t12 * t20;
	t5 = [0, -t19 * t3 + t26 * t4, (t8 * t25 - t4 * t7) * r_i_i_C(1) + (-t7 * t25 - t4 * t8) * r_i_i_C(2) + (t10 * t22 - t15 * t4) * pkin(3), t3, 0, 0; 0, -t19 * t1 + t26 * t2, (-t2 * t7 - t8 * t24) * r_i_i_C(1) + (-t2 * t8 + t7 * t24) * r_i_i_C(2) + (-t12 * t22 - t15 * t2) * pkin(3), t1, 0, 0; 1, (t26 * t16 + t19 * t18) * t11, (t13 * t8 - t7 * t23) * r_i_i_C(1) + (-t13 * t7 - t8 * t23) * r_i_i_C(2) + (t13 * t17 - t15 * t23) * pkin(3), -t11 * t18, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (125->30), mult. (200->52), div. (0->0), fcn. (249->10), ass. (0->29)
	t15 = qJ(3) + pkin(11);
	t13 = sin(t15);
	t14 = cos(t15);
	t23 = cos(qJ(3));
	t26 = r_i_i_C(3) + qJ(5);
	t34 = pkin(4) - r_i_i_C(2);
	t35 = pkin(3) * t23 + t26 * t13 + t34 * t14 + pkin(2);
	t33 = r_i_i_C(1) + qJ(4) + pkin(8);
	t16 = sin(pkin(10));
	t17 = sin(pkin(6));
	t32 = t16 * t17;
	t18 = cos(pkin(10));
	t31 = t17 * t18;
	t22 = sin(qJ(2));
	t30 = t17 * t22;
	t29 = t17 * t23;
	t19 = cos(pkin(6));
	t28 = t19 * t22;
	t24 = cos(qJ(2));
	t27 = t19 * t24;
	t21 = sin(qJ(3));
	t10 = -t16 * t28 + t18 * t24;
	t9 = t16 * t27 + t18 * t22;
	t8 = t16 * t24 + t18 * t28;
	t7 = t16 * t22 - t18 * t27;
	t5 = t13 * t30 - t14 * t19;
	t3 = t10 * t13 - t14 * t32;
	t1 = t13 * t8 + t14 * t31;
	t2 = [0, t33 * t10 - t35 * t9, t26 * (t10 * t14 + t13 * t32) - t34 * t3 + (-t10 * t21 + t16 * t29) * pkin(3), t9, t3, 0; 0, t33 * t8 - t35 * t7, t26 * (-t13 * t31 + t14 * t8) - t34 * t1 + (-t18 * t29 - t21 * t8) * pkin(3), t7, t1, 0; 1, (t33 * t22 + t35 * t24) * t17, t26 * (t13 * t19 + t14 * t30) - t34 * t5 + (t19 * t23 - t21 * t30) * pkin(3), -t17 * t24, t5, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (214->43), mult. (356->73), div. (0->0), fcn. (449->12), ass. (0->31)
	t17 = qJ(3) + pkin(11);
	t15 = sin(t17);
	t16 = cos(t17);
	t25 = cos(qJ(3));
	t21 = sin(qJ(6));
	t24 = cos(qJ(6));
	t29 = t21 * r_i_i_C(1) + t24 * r_i_i_C(2) + qJ(5);
	t33 = pkin(4) + pkin(9) + r_i_i_C(3);
	t39 = t25 * pkin(3) + t29 * t15 + t33 * t16 + pkin(2);
	t18 = sin(pkin(10));
	t19 = sin(pkin(6));
	t38 = t18 * t19;
	t23 = sin(qJ(2));
	t37 = t19 * t23;
	t26 = cos(qJ(2));
	t36 = t19 * t26;
	t35 = cos(pkin(6));
	t34 = cos(pkin(10));
	t32 = t18 * t35;
	t31 = t19 * t34;
	t30 = t35 * t34;
	t28 = t24 * r_i_i_C(1) - t21 * r_i_i_C(2) + pkin(5) + pkin(8) + qJ(4);
	t22 = sin(qJ(3));
	t10 = -t23 * t32 + t34 * t26;
	t9 = t34 * t23 + t26 * t32;
	t8 = t18 * t26 + t23 * t30;
	t7 = t18 * t23 - t26 * t30;
	t5 = t15 * t37 - t35 * t16;
	t3 = t10 * t15 - t16 * t38;
	t1 = t8 * t15 + t16 * t31;
	t2 = [0, t28 * t10 - t39 * t9, (-t10 * t22 + t25 * t38) * pkin(3) + t29 * (t10 * t16 + t15 * t38) - t33 * t3, t9, t3, (-t9 * t21 + t3 * t24) * r_i_i_C(1) + (-t3 * t21 - t9 * t24) * r_i_i_C(2); 0, t28 * t8 - t39 * t7, (-t22 * t8 - t25 * t31) * pkin(3) + t29 * (-t15 * t31 + t8 * t16) - t33 * t1, t7, t1, (t1 * t24 - t7 * t21) * r_i_i_C(1) + (-t1 * t21 - t7 * t24) * r_i_i_C(2); 1, (t28 * t23 + t39 * t26) * t19, (-t22 * t37 + t35 * t25) * pkin(3) + t29 * (t35 * t15 + t16 * t37) - t33 * t5, -t36, t5, (t21 * t36 + t5 * t24) * r_i_i_C(1) + (-t5 * t21 + t24 * t36) * r_i_i_C(2);];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end