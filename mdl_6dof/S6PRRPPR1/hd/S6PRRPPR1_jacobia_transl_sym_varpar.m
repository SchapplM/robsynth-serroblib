% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.21s
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
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (152->33), mult. (253->56), div. (0->0), fcn. (317->12), ass. (0->31)
	t17 = qJ(3) + pkin(11);
	t15 = sin(t17);
	t16 = cos(t17);
	t26 = cos(qJ(3));
	t18 = sin(pkin(12));
	t21 = cos(pkin(12));
	t30 = r_i_i_C(1) * t21 - r_i_i_C(2) * t18 + pkin(4);
	t34 = r_i_i_C(3) + qJ(5);
	t39 = t26 * pkin(3) + t34 * t15 + t30 * t16 + pkin(2);
	t19 = sin(pkin(10));
	t20 = sin(pkin(6));
	t38 = t19 * t20;
	t22 = cos(pkin(10));
	t37 = t20 * t22;
	t25 = sin(qJ(2));
	t36 = t20 * t25;
	t35 = t20 * t26;
	t33 = cos(pkin(6));
	t32 = t25 * t33;
	t27 = cos(qJ(2));
	t31 = t27 * t33;
	t29 = t18 * r_i_i_C(1) + t21 * r_i_i_C(2) + pkin(8) + qJ(4);
	t24 = sin(qJ(3));
	t10 = -t19 * t32 + t22 * t27;
	t9 = t19 * t31 + t22 * t25;
	t8 = t19 * t27 + t22 * t32;
	t7 = t19 * t25 - t22 * t31;
	t5 = t15 * t36 - t33 * t16;
	t3 = t10 * t15 - t16 * t38;
	t1 = t8 * t15 + t16 * t37;
	t2 = [0, t29 * t10 - t39 * t9, t34 * (t10 * t16 + t15 * t38) + (-t10 * t24 + t19 * t35) * pkin(3) - t30 * t3, t9, t3, 0; 0, t29 * t8 - t39 * t7, t34 * (-t15 * t37 + t8 * t16) + (-t22 * t35 - t24 * t8) * pkin(3) - t30 * t1, t7, t1, 0; 1, (t29 * t25 + t39 * t27) * t20, t34 * (t33 * t15 + t16 * t36) + (-t24 * t36 + t33 * t26) * pkin(3) - t30 * t5, -t20 * t27, t5, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (234->45), mult. (334->76), div. (0->0), fcn. (420->14), ass. (0->36)
	t21 = qJ(3) + pkin(11);
	t17 = sin(t21);
	t19 = cos(t21);
	t30 = cos(qJ(3));
	t20 = pkin(12) + qJ(6);
	t16 = sin(t20);
	t18 = cos(t20);
	t34 = t18 * r_i_i_C(1) - t16 * r_i_i_C(2) + cos(pkin(12)) * pkin(5) + pkin(4);
	t43 = r_i_i_C(3) + pkin(9) + qJ(5);
	t44 = t30 * pkin(3) + t43 * t17 + t34 * t19 + pkin(2);
	t23 = sin(pkin(10));
	t24 = sin(pkin(6));
	t42 = t23 * t24;
	t25 = cos(pkin(10));
	t41 = t24 * t25;
	t29 = sin(qJ(2));
	t40 = t24 * t29;
	t39 = t24 * t30;
	t31 = cos(qJ(2));
	t38 = t24 * t31;
	t37 = cos(pkin(6));
	t36 = t29 * t37;
	t35 = t31 * t37;
	t33 = sin(pkin(12)) * pkin(5) + t16 * r_i_i_C(1) + t18 * r_i_i_C(2) + qJ(4) + pkin(8);
	t28 = sin(qJ(3));
	t10 = -t23 * t36 + t25 * t31;
	t9 = t23 * t35 + t25 * t29;
	t8 = t23 * t31 + t25 * t36;
	t7 = t23 * t29 - t25 * t35;
	t6 = t37 * t17 + t19 * t40;
	t5 = t17 * t40 - t37 * t19;
	t4 = t10 * t19 + t17 * t42;
	t3 = t10 * t17 - t19 * t42;
	t2 = -t17 * t41 + t8 * t19;
	t1 = t8 * t17 + t19 * t41;
	t11 = [0, t33 * t10 - t44 * t9, t43 * t4 + (-t10 * t28 + t23 * t39) * pkin(3) - t34 * t3, t9, t3, (-t4 * t16 + t9 * t18) * r_i_i_C(1) + (-t9 * t16 - t4 * t18) * r_i_i_C(2); 0, t33 * t8 - t44 * t7, t43 * t2 + (-t25 * t39 - t28 * t8) * pkin(3) - t34 * t1, t7, t1, (-t2 * t16 + t7 * t18) * r_i_i_C(1) + (-t7 * t16 - t2 * t18) * r_i_i_C(2); 1, (t33 * t29 + t44 * t31) * t24, t43 * t6 + (-t28 * t40 + t37 * t30) * pkin(3) - t34 * t5, -t38, t5, (-t6 * t16 - t18 * t38) * r_i_i_C(1) + (t16 * t38 - t6 * t18) * r_i_i_C(2);];
	Ja_transl = t11;
else
	Ja_transl=NaN(3,6);
end