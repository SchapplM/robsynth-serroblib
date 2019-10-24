% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRR6
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:27
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRPRR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:35
	% EndTime: 2019-10-24 10:27:35
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:35
	% EndTime: 2019-10-24 10:27:35
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:35
	% EndTime: 2019-10-24 10:27:35
	% DurationCPUTime: 0.06s
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
	% StartTime: 2019-10-24 10:27:35
	% EndTime: 2019-10-24 10:27:35
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (25->11), mult. (63->20), div. (0->0), fcn. (78->8), ass. (0->13)
	t11 = cos(pkin(5));
	t12 = sin(qJ(2));
	t17 = t11 * t12;
	t13 = cos(qJ(2));
	t16 = t11 * t13;
	t15 = r_i_i_C(3) + qJ(3);
	t14 = r_i_i_C(1) * cos(pkin(10)) - r_i_i_C(2) * sin(pkin(10)) + pkin(2);
	t10 = cos(pkin(9));
	t8 = sin(pkin(5));
	t7 = sin(pkin(9));
	t3 = t10 * t12 + t7 * t16;
	t1 = -t10 * t16 + t7 * t12;
	t2 = [0, t15 * (t10 * t13 - t17 * t7) - t14 * t3, t3, 0, 0; 0, t15 * (t10 * t17 + t7 * t13) - t14 * t1, t1, 0, 0; 1, (t12 * t15 + t13 * t14) * t8, -t8 * t13, 0, 0;];
	Ja_transl = t2;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:35
	% EndTime: 2019-10-24 10:27:35
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (62->22), mult. (102->42), div. (0->0), fcn. (127->9), ass. (0->21)
	t23 = r_i_i_C(3) + pkin(7) + qJ(3);
	t10 = sin(pkin(9));
	t11 = sin(pkin(5));
	t22 = t10 * t11;
	t12 = cos(pkin(9));
	t21 = t11 * t12;
	t15 = sin(qJ(2));
	t20 = t11 * t15;
	t13 = cos(pkin(5));
	t19 = t13 * t15;
	t16 = cos(qJ(2));
	t18 = t13 * t16;
	t9 = pkin(10) + qJ(4);
	t7 = sin(t9);
	t8 = cos(t9);
	t17 = r_i_i_C(1) * t8 - r_i_i_C(2) * t7 + cos(pkin(10)) * pkin(3) + pkin(2);
	t4 = -t10 * t19 + t12 * t16;
	t3 = t10 * t18 + t12 * t15;
	t2 = t10 * t16 + t12 * t19;
	t1 = t10 * t15 - t12 * t18;
	t5 = [0, -t17 * t3 + t23 * t4, t3, (t8 * t22 - t4 * t7) * r_i_i_C(1) + (-t7 * t22 - t4 * t8) * r_i_i_C(2), 0; 0, -t17 * t1 + t23 * t2, t1, (-t2 * t7 - t8 * t21) * r_i_i_C(1) + (-t2 * t8 + t7 * t21) * r_i_i_C(2), 0; 1, (t23 * t15 + t17 * t16) * t11, -t11 * t16, (t13 * t8 - t7 * t20) * r_i_i_C(1) + (-t13 * t7 - t8 * t20) * r_i_i_C(2), 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:27:35
	% EndTime: 2019-10-24 10:27:35
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (168->36), mult. (279->64), div. (0->0), fcn. (353->11), ass. (0->29)
	t15 = pkin(10) + qJ(4);
	t13 = sin(t15);
	t14 = cos(t15);
	t19 = sin(qJ(5));
	t21 = cos(qJ(5));
	t25 = t21 * r_i_i_C(1) - t19 * r_i_i_C(2) + pkin(4);
	t34 = pkin(8) + r_i_i_C(3);
	t35 = t34 * t13 + t25 * t14 + cos(pkin(10)) * pkin(3) + pkin(2);
	t16 = sin(pkin(9));
	t17 = sin(pkin(5));
	t33 = t16 * t17;
	t20 = sin(qJ(2));
	t32 = t17 * t20;
	t22 = cos(qJ(2));
	t31 = t17 * t22;
	t30 = cos(pkin(5));
	t29 = cos(pkin(9));
	t28 = t16 * t30;
	t27 = t17 * t29;
	t26 = t30 * t29;
	t24 = t19 * r_i_i_C(1) + t21 * r_i_i_C(2) + pkin(7) + qJ(3);
	t10 = -t20 * t28 + t29 * t22;
	t9 = t29 * t20 + t22 * t28;
	t8 = t16 * t22 + t20 * t26;
	t7 = t16 * t20 - t22 * t26;
	t6 = t30 * t13 + t14 * t32;
	t4 = t10 * t14 + t13 * t33;
	t2 = -t13 * t27 + t8 * t14;
	t1 = [0, t24 * t10 - t35 * t9, t9, t34 * t4 + t25 * (-t10 * t13 + t14 * t33), (-t4 * t19 + t9 * t21) * r_i_i_C(1) + (-t9 * t19 - t4 * t21) * r_i_i_C(2); 0, t24 * t8 - t35 * t7, t7, t34 * t2 + t25 * (-t8 * t13 - t14 * t27), (-t2 * t19 + t7 * t21) * r_i_i_C(1) + (-t7 * t19 - t2 * t21) * r_i_i_C(2); 1, (t24 * t20 + t35 * t22) * t17, -t31, t34 * t6 + t25 * (-t13 * t32 + t30 * t14), (-t6 * t19 - t21 * t31) * r_i_i_C(1) + (t19 * t31 - t6 * t21) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end