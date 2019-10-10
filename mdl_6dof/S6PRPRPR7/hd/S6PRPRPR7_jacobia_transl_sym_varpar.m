% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:40
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
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
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.14s
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
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (74->23), mult. (186->40), div. (0->0), fcn. (235->8), ass. (0->23)
	t28 = pkin(4) - r_i_i_C(2);
	t13 = sin(pkin(6));
	t16 = sin(qJ(4));
	t27 = t13 * t16;
	t18 = cos(qJ(4));
	t26 = t13 * t18;
	t19 = cos(qJ(2));
	t25 = t13 * t19;
	t15 = cos(pkin(6));
	t17 = sin(qJ(2));
	t24 = t15 * t17;
	t23 = t15 * t19;
	t22 = r_i_i_C(3) + qJ(5);
	t21 = pkin(2) + pkin(8) + r_i_i_C(1);
	t20 = t28 * t16 - t22 * t18 + qJ(3);
	t14 = cos(pkin(10));
	t12 = sin(pkin(10));
	t9 = t15 * t16 + t18 * t25;
	t7 = t12 * t23 + t14 * t17;
	t5 = t12 * t17 - t14 * t23;
	t3 = t14 * t27 + t18 * t5;
	t1 = t12 * t27 - t18 * t7;
	t2 = [0, -t21 * t7 + t20 * (-t12 * t24 + t14 * t19), t7, t22 * (t12 * t26 + t16 * t7) - t28 * t1, t1, 0; 0, -t21 * t5 + t20 * (t12 * t19 + t14 * t24), t5, -t22 * (t14 * t26 - t16 * t5) + t28 * t3, -t3, 0; 1, (t20 * t17 + t21 * t19) * t13, -t25, -t28 * t9 + t22 * (t15 * t18 - t16 * t25), t9, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:40:56
	% EndTime: 2019-10-09 21:40:56
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (133->36), mult. (342->63), div. (0->0), fcn. (435->10), ass. (0->28)
	t15 = sin(pkin(6));
	t19 = sin(qJ(4));
	t33 = t15 * t19;
	t20 = sin(qJ(2));
	t32 = t15 * t20;
	t22 = cos(qJ(4));
	t31 = t15 * t22;
	t23 = cos(qJ(2));
	t30 = t15 * t23;
	t17 = cos(pkin(6));
	t29 = t17 * t20;
	t28 = t17 * t23;
	t27 = pkin(4) + pkin(9) + r_i_i_C(3);
	t18 = sin(qJ(6));
	t21 = cos(qJ(6));
	t26 = t18 * r_i_i_C(1) + t21 * r_i_i_C(2) + qJ(5);
	t25 = t21 * r_i_i_C(1) - t18 * r_i_i_C(2) + pkin(2) + pkin(5) + pkin(8);
	t24 = t27 * t19 - t26 * t22 + qJ(3);
	t16 = cos(pkin(10));
	t14 = sin(pkin(10));
	t11 = t17 * t19 + t22 * t30;
	t10 = -t14 * t29 + t16 * t23;
	t9 = t14 * t28 + t16 * t20;
	t8 = t14 * t23 + t16 * t29;
	t7 = t14 * t20 - t16 * t28;
	t3 = t16 * t33 + t7 * t22;
	t1 = t14 * t33 - t9 * t22;
	t2 = [0, t24 * t10 - t25 * t9, t9, t26 * (t14 * t31 + t9 * t19) - t27 * t1, t1, (t1 * t21 - t10 * t18) * r_i_i_C(1) + (-t1 * t18 - t10 * t21) * r_i_i_C(2); 0, t24 * t8 - t25 * t7, t7, -t26 * (t16 * t31 - t7 * t19) + t27 * t3, -t3, (-t8 * t18 - t3 * t21) * r_i_i_C(1) + (t3 * t18 - t8 * t21) * r_i_i_C(2); 1, (t24 * t20 + t25 * t23) * t15, -t30, t26 * (t17 * t22 - t19 * t30) - t27 * t11, t11, (t11 * t21 - t18 * t32) * r_i_i_C(1) + (-t11 * t18 - t21 * t32) * r_i_i_C(2);];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end