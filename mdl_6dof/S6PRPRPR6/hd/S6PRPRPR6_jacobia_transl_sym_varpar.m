% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
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
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (20->10), mult. (47->17), div. (0->0), fcn. (60->6), ass. (0->14)
	t8 = cos(pkin(6));
	t9 = sin(qJ(2));
	t15 = t8 * t9;
	t14 = pkin(2) - r_i_i_C(2);
	t10 = cos(qJ(2));
	t5 = sin(pkin(10));
	t13 = t5 * t10;
	t7 = cos(pkin(10));
	t12 = t7 * t10;
	t11 = r_i_i_C(3) + qJ(3);
	t6 = sin(pkin(6));
	t3 = t13 * t8 + t7 * t9;
	t1 = -t12 * t8 + t5 * t9;
	t2 = [0, t11 * (-t15 * t5 + t12) - t14 * t3, t3, 0, 0, 0; 0, t11 * (t15 * t7 + t13) - t14 * t1, t1, 0, 0, 0; 1, (t10 * t14 + t11 * t9) * t6, -t6 * t10, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (43->20), mult. (109->39), div. (0->0), fcn. (136->8), ass. (0->19)
	t10 = sin(qJ(4));
	t7 = sin(pkin(6));
	t21 = t10 * t7;
	t11 = sin(qJ(2));
	t6 = sin(pkin(10));
	t20 = t11 * t6;
	t12 = cos(qJ(4));
	t19 = t12 * t7;
	t13 = cos(qJ(2));
	t9 = cos(pkin(6));
	t18 = t13 * t9;
	t17 = t7 * t13;
	t8 = cos(pkin(10));
	t16 = t8 * t11;
	t15 = pkin(2) + pkin(8) + r_i_i_C(3);
	t14 = r_i_i_C(1) * t10 + r_i_i_C(2) * t12 + qJ(3);
	t3 = t6 * t18 + t16;
	t1 = -t8 * t18 + t20;
	t2 = [0, t14 * (t13 * t8 - t9 * t20) - t15 * t3, t3, (t12 * t3 - t6 * t21) * r_i_i_C(1) + (-t10 * t3 - t6 * t19) * r_i_i_C(2), 0, 0; 0, t14 * (t13 * t6 + t9 * t16) - t15 * t1, t1, (t1 * t12 + t8 * t21) * r_i_i_C(1) + (-t1 * t10 + t8 * t19) * r_i_i_C(2), 0, 0; 1, (t14 * t11 + t15 * t13) * t7, -t17, (-t9 * t10 - t12 * t17) * r_i_i_C(1) + (t10 * t17 - t9 * t12) * r_i_i_C(2), 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (92->26), mult. (239->44), div. (0->0), fcn. (303->10), ass. (0->25)
	t15 = sin(pkin(6));
	t19 = sin(qJ(4));
	t31 = t15 * t19;
	t21 = cos(qJ(4));
	t30 = t15 * t21;
	t22 = cos(qJ(2));
	t29 = t15 * t22;
	t18 = cos(pkin(6));
	t20 = sin(qJ(2));
	t28 = t18 * t20;
	t27 = t18 * t22;
	t26 = r_i_i_C(3) + qJ(5);
	t13 = sin(pkin(11));
	t16 = cos(pkin(11));
	t25 = r_i_i_C(1) * t16 - r_i_i_C(2) * t13 + pkin(4);
	t24 = t13 * r_i_i_C(1) + t16 * r_i_i_C(2) + pkin(2) + pkin(8);
	t23 = t25 * t19 - t26 * t21 + qJ(3);
	t17 = cos(pkin(10));
	t14 = sin(pkin(10));
	t10 = t18 * t19 + t21 * t29;
	t8 = t14 * t27 + t17 * t20;
	t6 = t14 * t20 - t17 * t27;
	t3 = t17 * t31 + t6 * t21;
	t1 = t14 * t31 - t8 * t21;
	t2 = [0, -t24 * t8 + t23 * (-t14 * t28 + t17 * t22), t8, t26 * (t14 * t30 + t8 * t19) - t25 * t1, t1, 0; 0, -t24 * t6 + t23 * (t14 * t22 + t17 * t28), t6, -t26 * (t17 * t30 - t6 * t19) + t25 * t3, -t3, 0; 1, (t23 * t20 + t24 * t22) * t15, -t29, t26 * (t18 * t21 - t19 * t29) - t25 * t10, t10, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (162->38), mult. (320->65), div. (0->0), fcn. (406->12), ass. (0->32)
	t37 = r_i_i_C(3) + pkin(9) + qJ(5);
	t20 = sin(pkin(6));
	t24 = sin(qJ(4));
	t36 = t20 * t24;
	t25 = sin(qJ(2));
	t35 = t20 * t25;
	t26 = cos(qJ(4));
	t34 = t20 * t26;
	t27 = cos(qJ(2));
	t33 = t20 * t27;
	t22 = cos(pkin(6));
	t32 = t22 * t25;
	t31 = t22 * t27;
	t17 = pkin(11) + qJ(6);
	t15 = sin(t17);
	t16 = cos(t17);
	t30 = t16 * r_i_i_C(1) - t15 * r_i_i_C(2) + cos(pkin(11)) * pkin(5) + pkin(4);
	t29 = sin(pkin(11)) * pkin(5) + t15 * r_i_i_C(1) + t16 * r_i_i_C(2) + pkin(2) + pkin(8);
	t28 = t30 * t24 - t37 * t26 + qJ(3);
	t21 = cos(pkin(10));
	t19 = sin(pkin(10));
	t12 = t22 * t26 - t24 * t33;
	t11 = t22 * t24 + t26 * t33;
	t10 = -t19 * t32 + t21 * t27;
	t9 = t19 * t31 + t21 * t25;
	t8 = t19 * t27 + t21 * t32;
	t7 = t19 * t25 - t21 * t31;
	t4 = t21 * t34 - t7 * t24;
	t3 = t21 * t36 + t7 * t26;
	t2 = t19 * t34 + t9 * t24;
	t1 = t19 * t36 - t9 * t26;
	t5 = [0, t28 * t10 - t29 * t9, t9, -t30 * t1 + t37 * t2, t1, (t10 * t16 - t2 * t15) * r_i_i_C(1) + (-t10 * t15 - t2 * t16) * r_i_i_C(2); 0, t28 * t8 - t29 * t7, t7, t30 * t3 - t37 * t4, -t3, (t4 * t15 + t8 * t16) * r_i_i_C(1) + (-t8 * t15 + t4 * t16) * r_i_i_C(2); 1, (t28 * t25 + t29 * t27) * t20, -t33, -t30 * t11 + t37 * t12, t11, (-t12 * t15 + t16 * t35) * r_i_i_C(1) + (-t12 * t16 - t15 * t35) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end