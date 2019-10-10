% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:56
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
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
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->5), mult. (8->4), div. (0->0), fcn. (10->2), ass. (0->5)
	t4 = pkin(1) - r_i_i_C(2);
	t3 = r_i_i_C(3) + qJ(2);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t1 * t4 + t2 * t3, t1, 0, 0, 0, 0; t1 * t3 + t2 * t4, -t2, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (9->5), mult. (10->4), div. (0->0), fcn. (14->2), ass. (0->5)
	t4 = r_i_i_C(2) + qJ(2);
	t3 = pkin(1) + r_i_i_C(3) + qJ(3);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t3 * t1 + t4 * t2, t1, t2, 0, 0, 0; t4 * t1 + t3 * t2, -t2, t1, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (16->7), mult. (28->10), div. (0->0), fcn. (32->4), ass. (0->9)
	t8 = -pkin(7) - r_i_i_C(3) + qJ(2);
	t1 = sin(qJ(4));
	t3 = cos(qJ(4));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t5 = pkin(1) + qJ(3) - t6;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t5 * t2 + t8 * t4, t2, t4, t7 * t4, 0, 0; t8 * t2 + t5 * t4, -t4, t2, t7 * t2, 0, 0; 0, 0, 0, t6, 0, 0;];
	Ja_transl = t9;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (40->23), mult. (87->32), div. (0->0), fcn. (99->6), ass. (0->21)
	t18 = pkin(8) + r_i_i_C(3);
	t9 = cos(qJ(4));
	t14 = t18 * t9;
	t6 = sin(qJ(4));
	t21 = -t6 * pkin(4) - pkin(1) - qJ(3) + t14;
	t5 = sin(qJ(5));
	t7 = sin(qJ(1));
	t20 = t7 * t5;
	t8 = cos(qJ(5));
	t19 = t7 * t8;
	t10 = cos(qJ(1));
	t17 = t10 * t5;
	t16 = t10 * t8;
	t15 = -pkin(7) + qJ(2);
	t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5 + pkin(4);
	t11 = t12 * t9 + t18 * t6;
	t4 = t6 * t16 - t20;
	t3 = -t6 * t17 - t19;
	t2 = -t6 * t19 - t17;
	t1 = t6 * t20 - t16;
	t13 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * t10 + t21 * t7, t7, t10, t11 * t10, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t21 * t10 + t15 * t7, -t10, t7, t11 * t7, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, 0, 0, -t12 * t6 + t14, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t9, 0;];
	Ja_transl = t13;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:56:17
	% EndTime: 2019-10-09 23:56:17
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (66->25), mult. (146->33), div. (0->0), fcn. (171->6), ass. (0->23)
	t10 = cos(qJ(4));
	t19 = pkin(8) + r_i_i_C(2);
	t14 = t19 * t10;
	t7 = sin(qJ(4));
	t24 = -pkin(4) * t7 - pkin(1) - qJ(3) + t14;
	t15 = r_i_i_C(3) + qJ(6);
	t20 = pkin(5) + r_i_i_C(1);
	t6 = sin(qJ(5));
	t9 = cos(qJ(5));
	t23 = t15 * t6 + t20 * t9 + pkin(4);
	t8 = sin(qJ(1));
	t22 = t8 * t6;
	t21 = t8 * t9;
	t11 = cos(qJ(1));
	t18 = t11 * t6;
	t17 = t11 * t9;
	t16 = -pkin(7) + qJ(2);
	t12 = t23 * t10 + t19 * t7;
	t4 = t7 * t17 - t22;
	t3 = t7 * t18 + t21;
	t2 = t7 * t21 + t18;
	t1 = t7 * t22 - t17;
	t5 = [-t15 * t1 + t16 * t11 - t20 * t2 + t24 * t8, t8, t11, t12 * t11, t15 * t4 - t20 * t3, t3; -t24 * t11 + t15 * t3 + t16 * t8 + t20 * t4, -t11, t8, t12 * t8, -t20 * t1 + t15 * t2, t1; 0, 0, 0, -t23 * t7 + t14, (t15 * t9 - t20 * t6) * t10, t10 * t6;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end