% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:26
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPPRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:18
	% EndTime: 2019-10-09 21:26:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:18
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:18
	% EndTime: 2019-10-09 21:26:19
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
	% StartTime: 2019-10-09 21:26:18
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (23->15), mult. (60->31), div. (0->0), fcn. (77->8), ass. (0->14)
	t10 = cos(pkin(6));
	t12 = cos(qJ(2));
	t13 = t10 * t12;
	t11 = sin(qJ(2));
	t5 = sin(pkin(11));
	t8 = cos(pkin(11));
	t4 = -t11 * t8 - t12 * t5;
	t3 = t11 * t5 - t12 * t8;
	t9 = cos(pkin(10));
	t7 = sin(pkin(6));
	t6 = sin(pkin(10));
	t2 = t4 * t10;
	t1 = t3 * t10;
	t14 = [0, (t6 * t1 + t9 * t4) * r_i_i_C(1) + (-t6 * t2 + t9 * t3) * r_i_i_C(2) + (-t9 * t11 - t6 * t13) * pkin(2), t6 * t7, 0, 0, 0; 0, (-t9 * t1 + t6 * t4) * r_i_i_C(1) + (t9 * t2 + t6 * t3) * r_i_i_C(2) + (-t6 * t11 + t9 * t13) * pkin(2), -t9 * t7, 0, 0, 0; 1, (t12 * pkin(2) - t3 * r_i_i_C(1) + t4 * r_i_i_C(2)) * t7, t10, 0, 0, 0;];
	Ja_transl = t14;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:18
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (50->19), mult. (126->32), div. (0->0), fcn. (167->8), ass. (0->19)
	t23 = pkin(3) - r_i_i_C(2);
	t16 = cos(pkin(6));
	t18 = cos(qJ(2));
	t22 = t16 * t18;
	t21 = r_i_i_C(3) + qJ(4);
	t11 = sin(pkin(11));
	t14 = cos(pkin(11));
	t17 = sin(qJ(2));
	t20 = t11 * t18 + t17 * t14;
	t9 = t17 * t11 - t14 * t18;
	t19 = t9 * t16;
	t15 = cos(pkin(10));
	t13 = sin(pkin(6));
	t12 = sin(pkin(10));
	t8 = t20 * t16;
	t6 = t9 * t13;
	t4 = t12 * t19 - t15 * t20;
	t2 = -t12 * t20 - t15 * t19;
	t1 = [0, t23 * t4 - t21 * (t12 * t8 + t15 * t9) + (-t12 * t22 - t15 * t17) * pkin(2), t12 * t13, -t4, 0, 0; 0, t23 * t2 - t21 * (t12 * t9 - t15 * t8) + (-t12 * t17 + t15 * t22) * pkin(2), -t15 * t13, -t2, 0, 0; 1, -t23 * t6 + (pkin(2) * t18 + t20 * t21) * t13, t16, t6, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:18
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (93->29), mult. (240->52), div. (0->0), fcn. (315->10), ass. (0->23)
	t12 = sin(pkin(11));
	t15 = cos(pkin(11));
	t19 = sin(qJ(2));
	t21 = cos(qJ(2));
	t9 = t19 * t12 - t21 * t15;
	t13 = sin(pkin(10));
	t14 = sin(pkin(6));
	t31 = t13 * t14;
	t16 = cos(pkin(10));
	t30 = t16 * t14;
	t17 = cos(pkin(6));
	t29 = t17 * t21;
	t26 = pkin(3) + pkin(8) + r_i_i_C(3);
	t25 = t21 * t12 + t19 * t15;
	t18 = sin(qJ(5));
	t20 = cos(qJ(5));
	t24 = t18 * r_i_i_C(1) + t20 * r_i_i_C(2) + qJ(4);
	t23 = t25 * t17;
	t22 = t9 * t17;
	t7 = t9 * t14;
	t4 = t13 * t22 - t16 * t25;
	t2 = -t13 * t25 - t16 * t22;
	t1 = [0, (-t13 * t29 - t16 * t19) * pkin(2) + t26 * t4 - t24 * (t13 * t23 + t16 * t9), t31, -t4, (-t18 * t31 - t4 * t20) * r_i_i_C(1) + (t4 * t18 - t20 * t31) * r_i_i_C(2), 0; 0, (-t13 * t19 + t16 * t29) * pkin(2) + t26 * t2 - t24 * (t13 * t9 - t16 * t23), -t30, -t2, (t18 * t30 - t2 * t20) * r_i_i_C(1) + (t2 * t18 + t20 * t30) * r_i_i_C(2), 0; 1, -t26 * t7 + (pkin(2) * t21 + t24 * t25) * t14, t17, t7, (-t17 * t18 + t7 * t20) * r_i_i_C(1) + (-t17 * t20 - t7 * t18) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:18
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (215->45), mult. (560->78), div. (0->0), fcn. (739->12), ass. (0->35)
	t33 = cos(qJ(2));
	t41 = cos(pkin(11));
	t24 = sin(pkin(11));
	t30 = sin(qJ(2));
	t44 = t30 * t24;
	t18 = -t33 * t41 + t44;
	t25 = sin(pkin(10));
	t27 = cos(pkin(10));
	t42 = cos(pkin(6));
	t38 = t42 * t41;
	t40 = t33 * t42;
	t43 = -t24 * t40 - t30 * t38;
	t11 = -t27 * t18 + t25 * t43;
	t51 = t25 * t18 + t27 * t43;
	t29 = sin(qJ(5));
	t32 = cos(qJ(5));
	t28 = sin(qJ(6));
	t31 = cos(qJ(6));
	t37 = t31 * r_i_i_C(1) - t28 * r_i_i_C(2) + pkin(5);
	t49 = pkin(9) + r_i_i_C(3);
	t50 = t37 * t29 - t49 * t32 + qJ(4);
	t26 = sin(pkin(6));
	t47 = t25 * t26;
	t45 = t27 * t26;
	t36 = t28 * r_i_i_C(1) + t31 * r_i_i_C(2) + pkin(3) + pkin(8);
	t19 = -t33 * t24 - t30 * t41;
	t35 = t33 * t38 - t42 * t44;
	t17 = t19 * t26;
	t16 = t18 * t26;
	t13 = t16 * t29 + t42 * t32;
	t10 = t27 * t19 - t25 * t35;
	t7 = t25 * t19 + t27 * t35;
	t4 = t7 * t29 + t32 * t45;
	t2 = -t10 * t29 + t32 * t47;
	t1 = [0, (-t25 * t40 - t27 * t30) * pkin(2) + t36 * t10 + t50 * t11, t47, -t10, t49 * t2 + t37 * (-t10 * t32 - t29 * t47), (t11 * t31 - t2 * t28) * r_i_i_C(1) + (-t11 * t28 - t2 * t31) * r_i_i_C(2); 0, (-t25 * t30 + t27 * t40) * pkin(2) + t36 * t7 - t50 * t51, -t45, -t7, -t49 * t4 + t37 * (t29 * t45 - t7 * t32), (t4 * t28 - t31 * t51) * r_i_i_C(1) + (t28 * t51 + t4 * t31) * r_i_i_C(2); 1, t26 * t33 * pkin(2) - t36 * t16 - t50 * t17, t42, t16, t49 * t13 + t37 * (t16 * t32 - t42 * t29), (-t13 * t28 - t17 * t31) * r_i_i_C(1) + (-t13 * t31 + t17 * t28) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end