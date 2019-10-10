% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:23
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPPR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPPR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:23
	% EndTime: 2019-10-10 09:23:23
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:23
	% EndTime: 2019-10-10 09:23:24
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
	% StartTime: 2019-10-10 09:23:23
	% EndTime: 2019-10-10 09:23:24
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
	% StartTime: 2019-10-10 09:23:23
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (22->9), mult. (44->12), div. (0->0), fcn. (47->4), ass. (0->11)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(3) + qJ(3);
	t9 = pkin(2) - r_i_i_C(2);
	t6 = t7 * t1 + t9 * t3;
	t10 = pkin(1) + t6;
	t8 = pkin(7) + r_i_i_C(1);
	t5 = -t9 * t1 + t7 * t3;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t11 = [-t10 * t2 + t8 * t4, t5 * t4, t4 * t1, 0, 0, 0; t10 * t4 + t8 * t2, t5 * t2, t2 * t1, 0, 0, 0; 0, t6, -t3, 0, 0, 0;];
	Ja_transl = t11;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:23
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (36->12), mult. (78->18), div. (0->0), fcn. (88->6), ass. (0->13)
	t1 = sin(pkin(9));
	t2 = cos(pkin(9));
	t10 = r_i_i_C(1) * t1 + r_i_i_C(2) * t2 + qJ(3);
	t11 = pkin(2) + r_i_i_C(3) + qJ(4);
	t3 = sin(qJ(2));
	t5 = cos(qJ(2));
	t8 = t10 * t3 + t11 * t5;
	t12 = pkin(1) + t8;
	t9 = r_i_i_C(1) * t2 - r_i_i_C(2) * t1 + pkin(3) + pkin(7);
	t7 = t10 * t5 - t11 * t3;
	t6 = cos(qJ(1));
	t4 = sin(qJ(1));
	t13 = [-t12 * t4 + t9 * t6, t7 * t6, t6 * t3, t6 * t5, 0, 0; t12 * t6 + t9 * t4, t7 * t4, t4 * t3, t4 * t5, 0, 0; 0, t8, -t5, t3, 0, 0;];
	Ja_transl = t13;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:23
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (52->22), mult. (117->30), div. (0->0), fcn. (136->6), ass. (0->19)
	t15 = pkin(2) + r_i_i_C(2) + qJ(4);
	t9 = cos(qJ(2));
	t13 = t15 * t9;
	t7 = sin(qJ(2));
	t21 = t7 * qJ(3) + pkin(1) + t13;
	t8 = sin(qJ(1));
	t20 = t8 * t7;
	t19 = pkin(3) + pkin(7);
	t18 = pkin(4) + r_i_i_C(1);
	t10 = cos(qJ(1));
	t17 = t10 * t7;
	t16 = r_i_i_C(3) + qJ(5);
	t5 = sin(pkin(9));
	t6 = cos(pkin(9));
	t12 = -t16 * t6 + t18 * t5 + qJ(3);
	t11 = t12 * t9 - t15 * t7;
	t3 = t10 * t5 + t6 * t20;
	t1 = -t6 * t17 + t8 * t5;
	t2 = [t18 * (t10 * t6 - t5 * t20) + t16 * t3 + t19 * t10 - t21 * t8, t11 * t10, t17, t10 * t9, t1, 0; t19 * t8 + t18 * (t5 * t17 + t8 * t6) + t16 * t1 + t21 * t10, t11 * t8, t20, t8 * t9, -t3, 0; 0, t12 * t7 + t13, -t9, t7, t9 * t6, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:23:23
	% EndTime: 2019-10-10 09:23:24
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (95->36), mult. (230->55), div. (0->0), fcn. (279->8), ass. (0->26)
	t11 = sin(qJ(2));
	t14 = cos(qJ(2));
	t22 = r_i_i_C(3) + pkin(8) - qJ(4) - pkin(2);
	t20 = t22 * t14;
	t27 = -t11 * qJ(3) - pkin(1) + t20;
	t26 = pkin(3) + pkin(7);
	t25 = pkin(4) + pkin(5);
	t12 = sin(qJ(1));
	t24 = t12 * t11;
	t15 = cos(qJ(1));
	t23 = t15 * t11;
	t10 = sin(qJ(6));
	t13 = cos(qJ(6));
	t19 = t10 * r_i_i_C(1) + t13 * r_i_i_C(2) + qJ(5);
	t18 = t13 * r_i_i_C(1) - t10 * r_i_i_C(2) + t25;
	t8 = sin(pkin(9));
	t9 = cos(pkin(9));
	t17 = t18 * t8 - t19 * t9 + qJ(3);
	t16 = t22 * t11 + t17 * t14;
	t6 = t15 * t9 - t8 * t24;
	t5 = t15 * t8 + t9 * t24;
	t4 = t12 * t9 + t8 * t23;
	t3 = t12 * t8 - t9 * t23;
	t2 = t3 * t10 + t4 * t13;
	t1 = -t4 * t10 + t3 * t13;
	t7 = [t27 * t12 + t26 * t15 + t18 * t6 + t19 * t5, t16 * t15, t23, t15 * t14, t3, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t3 * qJ(5) + t26 * t12 - t27 * t15 + t25 * t4, t16 * t12, t24, t12 * t14, -t5, (t6 * t10 - t5 * t13) * r_i_i_C(1) + (t5 * t10 + t6 * t13) * r_i_i_C(2); 0, t17 * t11 - t20, -t14, t11, t14 * t9, ((t10 * t8 + t13 * t9) * r_i_i_C(1) + (-t10 * t9 + t13 * t8) * r_i_i_C(2)) * t14;];
	Ja_transl = t7;
else
	Ja_transl=NaN(3,6);
end