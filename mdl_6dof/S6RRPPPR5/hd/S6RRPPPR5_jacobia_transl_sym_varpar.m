% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:25
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPPR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t2 * t5 + t4 * t8, t6 * t4, 0, 0, 0, 0; t2 * t8 + t4 * t5, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (29->11), mult. (65->16), div. (0->0), fcn. (72->6), ass. (0->13)
	t1 = sin(pkin(9));
	t2 = cos(pkin(9));
	t10 = r_i_i_C(1) * t2 - r_i_i_C(2) * t1 + pkin(2);
	t11 = r_i_i_C(3) + qJ(3);
	t3 = sin(qJ(2));
	t5 = cos(qJ(2));
	t8 = t10 * t5 + t11 * t3;
	t12 = pkin(1) + t8;
	t9 = t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + pkin(7);
	t7 = -t10 * t3 + t11 * t5;
	t6 = cos(qJ(1));
	t4 = sin(qJ(1));
	t13 = [-t12 * t4 + t9 * t6, t7 * t6, t6 * t3, 0, 0, 0; t12 * t6 + t9 * t4, t7 * t4, t4 * t3, 0, 0, 0; 0, t8, -t5, 0, 0, 0;];
	Ja_transl = t13;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (45->20), mult. (104->30), div. (0->0), fcn. (120->6), ass. (0->18)
	t15 = r_i_i_C(1) + qJ(3);
	t7 = sin(qJ(2));
	t12 = t15 * t7;
	t9 = cos(qJ(2));
	t20 = t9 * pkin(2) + pkin(1) + t12;
	t14 = r_i_i_C(3) + qJ(4);
	t17 = pkin(3) - r_i_i_C(2);
	t5 = sin(pkin(9));
	t6 = cos(pkin(9));
	t19 = t14 * t5 + t17 * t6 + pkin(2);
	t8 = sin(qJ(1));
	t18 = t8 * t9;
	t10 = cos(qJ(1));
	t16 = t10 * t9;
	t11 = t15 * t9 - t19 * t7;
	t3 = t5 * t16 - t8 * t6;
	t1 = t10 * t6 + t5 * t18;
	t2 = [pkin(7) * t10 - t17 * (-t10 * t5 + t6 * t18) - t14 * t1 - t20 * t8, t11 * t10, t10 * t7, t3, 0, 0; t8 * pkin(7) + t17 * (t6 * t16 + t8 * t5) + t14 * t3 + t20 * t10, t11 * t8, t8 * t7, t1, 0, 0; 0, t19 * t9 + t12, -t9, t7 * t5, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (59->20), mult. (136->31), div. (0->0), fcn. (159->6), ass. (0->20)
	t10 = cos(qJ(2));
	t15 = pkin(4) - r_i_i_C(2) + qJ(3);
	t8 = sin(qJ(2));
	t13 = t15 * t8;
	t21 = t10 * pkin(2) + pkin(1) + t13;
	t16 = pkin(3) + r_i_i_C(3) + qJ(5);
	t17 = r_i_i_C(1) + qJ(4);
	t6 = sin(pkin(9));
	t7 = cos(pkin(9));
	t20 = t16 * t7 + t17 * t6 + pkin(2);
	t9 = sin(qJ(1));
	t19 = t10 * t9;
	t11 = cos(qJ(1));
	t18 = t10 * t11;
	t12 = t15 * t10 - t20 * t8;
	t4 = t7 * t18 + t9 * t6;
	t3 = t6 * t18 - t9 * t7;
	t2 = -t11 * t6 + t7 * t19;
	t1 = t11 * t7 + t6 * t19;
	t5 = [pkin(7) * t11 - t17 * t1 - t16 * t2 - t21 * t9, t12 * t11, t11 * t8, t3, t4, 0; t9 * pkin(7) + t21 * t11 + t16 * t4 + t17 * t3, t12 * t9, t9 * t8, t1, t2, 0; 0, t20 * t10 + t13, -t10, t8 * t6, t8 * t7, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:25:09
	% EndTime: 2019-10-10 09:25:09
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (102->35), mult. (249->56), div. (0->0), fcn. (302->8), ass. (0->26)
	t14 = cos(qJ(2));
	t11 = sin(qJ(2));
	t20 = r_i_i_C(3) + pkin(8) + pkin(4) + qJ(3);
	t19 = t20 * t11;
	t27 = t14 * pkin(2) + pkin(1) + t19;
	t10 = sin(qJ(6));
	t13 = cos(qJ(6));
	t22 = pkin(5) + qJ(4);
	t17 = t13 * r_i_i_C(1) - t10 * r_i_i_C(2) + t22;
	t23 = pkin(3) + qJ(5);
	t18 = t10 * r_i_i_C(1) + t13 * r_i_i_C(2) + t23;
	t8 = sin(pkin(9));
	t9 = cos(pkin(9));
	t26 = t17 * t8 + t18 * t9 + pkin(2);
	t12 = sin(qJ(1));
	t25 = t12 * t14;
	t15 = cos(qJ(1));
	t24 = t14 * t15;
	t16 = -t26 * t11 + t20 * t14;
	t6 = t12 * t8 + t9 * t24;
	t5 = -t12 * t9 + t8 * t24;
	t4 = -t15 * t8 + t9 * t25;
	t3 = t15 * t9 + t8 * t25;
	t2 = t6 * t10 + t5 * t13;
	t1 = -t5 * t10 + t6 * t13;
	t7 = [t15 * pkin(7) - t27 * t12 - t17 * t3 - t18 * t4, t16 * t15, t15 * t11, t5, t6, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t12 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t27 * t15 + t22 * t5 + t23 * t6, t16 * t12, t12 * t11, t3, t4, (-t3 * t10 + t4 * t13) * r_i_i_C(1) + (-t4 * t10 - t3 * t13) * r_i_i_C(2); 0, t26 * t14 + t19, -t14, t11 * t8, t11 * t9, ((-t10 * t8 + t13 * t9) * r_i_i_C(1) + (-t10 * t9 - t13 * t8) * r_i_i_C(2)) * t11;];
	Ja_transl = t7;
else
	Ja_transl=NaN(3,6);
end