% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR1
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:57
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobia_transl_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:57:44
	% EndTime: 2019-10-09 20:57:44
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:57:44
	% EndTime: 2019-10-09 20:57:44
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:57:44
	% EndTime: 2019-10-09 20:57:44
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (5->4), mult. (6->4), div. (0->0), fcn. (8->2), ass. (0->4)
	t3 = r_i_i_C(3) + qJ(2);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t4 = [-t1 * r_i_i_C(1) + t2 * t3, t1, 0, 0, 0; t2 * r_i_i_C(1) + t1 * t3, -t2, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:57:44
	% EndTime: 2019-10-09 20:57:44
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (10->6), mult. (22->10), div. (0->0), fcn. (24->4), ass. (0->8)
	t7 = r_i_i_C(3) + qJ(2);
	t1 = sin(qJ(3));
	t3 = cos(qJ(3));
	t6 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t5 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t8 = [-t6 * t2 + t7 * t4, t2, t5 * t4, 0, 0; t7 * t2 + t6 * t4, -t4, t5 * t2, 0, 0; 0, 0, t6, 0, 0;];
	Ja_transl = t8;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:57:44
	% EndTime: 2019-10-09 20:57:44
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (24->17), mult. (63->33), div. (0->0), fcn. (73->6), ass. (0->16)
	t6 = sin(qJ(3));
	t15 = t6 * r_i_i_C(3);
	t7 = sin(qJ(1));
	t9 = cos(qJ(3));
	t14 = t7 * t9;
	t10 = cos(qJ(1));
	t13 = t10 * t9;
	t5 = sin(qJ(4));
	t8 = cos(qJ(4));
	t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5;
	t11 = r_i_i_C(3) * t9 - t12 * t6;
	t4 = t8 * t13 + t7 * t5;
	t3 = -t5 * t13 + t7 * t8;
	t2 = t10 * t5 - t8 * t14;
	t1 = t10 * t8 + t5 * t14;
	t16 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t10 * qJ(2) - t7 * t15, t7, t11 * t10, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t7 * qJ(2) + t10 * t15, -t10, t11 * t7, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0; 0, 0, t12 * t9 + t15, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6, 0;];
	Ja_transl = t16;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:57:44
	% EndTime: 2019-10-09 20:57:44
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (54->30), mult. (150->59), div. (0->0), fcn. (184->8), ass. (0->30)
	t10 = sin(qJ(4));
	t31 = t10 * r_i_i_C(3);
	t11 = sin(qJ(3));
	t9 = sin(qJ(5));
	t30 = t11 * t9;
	t13 = cos(qJ(5));
	t29 = t11 * t13;
	t14 = cos(qJ(4));
	t28 = t11 * t14;
	t16 = cos(qJ(1));
	t27 = t11 * t16;
	t12 = sin(qJ(1));
	t26 = t12 * t10;
	t15 = cos(qJ(3));
	t25 = t14 * t15;
	t24 = t14 * t16;
	t23 = t16 * t10;
	t22 = t13 * r_i_i_C(1) - t9 * r_i_i_C(2);
	t4 = t12 * t25 - t23;
	t21 = -t12 * t30 - t13 * t4;
	t20 = t12 * t29 - t4 * t9;
	t19 = t13 * t15 + t9 * t28;
	t18 = -t13 * t28 + t15 * t9;
	t17 = t18 * r_i_i_C(1) + t19 * r_i_i_C(2) - t11 * t31;
	t6 = t15 * t24 + t26;
	t5 = -t12 * t14 + t15 * t23;
	t3 = -t15 * t26 - t24;
	t2 = t6 * t13 + t9 * t27;
	t1 = t13 * t27 - t6 * t9;
	t7 = [t21 * r_i_i_C(1) - t20 * r_i_i_C(2) + t3 * r_i_i_C(3) + t16 * qJ(2), t12, t17 * t16, r_i_i_C(3) * t6 - t22 * t5, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; r_i_i_C(1) * t2 + r_i_i_C(2) * t1 + r_i_i_C(3) * t5 + qJ(2) * t12, -t16, t17 * t12, r_i_i_C(3) * t4 + t22 * t3, t20 * r_i_i_C(1) + t21 * r_i_i_C(2); 0, 0, (t13 * t25 + t30) * r_i_i_C(1) + (-t9 * t25 + t29) * r_i_i_C(2) + t15 * t31, (r_i_i_C(3) * t14 - t22 * t10) * t11, -t19 * r_i_i_C(1) + t18 * r_i_i_C(2);];
	Ja_transl = t7;
else
	Ja_transl=NaN(3,5);
end