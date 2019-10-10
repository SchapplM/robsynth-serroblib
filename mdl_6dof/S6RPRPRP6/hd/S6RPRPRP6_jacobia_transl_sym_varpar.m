% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
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
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(9)) - r_i_i_C(2) * sin(pkin(9)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (26->9), mult. (26->11), div. (0->0), fcn. (28->5), ass. (0->10)
	t11 = r_i_i_C(3) + pkin(7) + qJ(2);
	t4 = pkin(9) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = cos(pkin(9)) * pkin(2) + pkin(1) + t10;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t8 * t6, t6, t9 * t7, 0, 0, 0; t11 * t6 + t8 * t7, -t7, t9 * t6, 0, 0, 0; 0, 0, t10, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (50->12), mult. (46->13), div. (0->0), fcn. (51->5), ass. (0->12)
	t10 = r_i_i_C(3) + qJ(4);
	t12 = pkin(3) - r_i_i_C(2);
	t4 = pkin(9) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t9 = t10 * t2 + t12 * t3;
	t13 = cos(pkin(9)) * pkin(2) + pkin(1) + t9;
	t11 = r_i_i_C(1) + pkin(7) + qJ(2);
	t8 = t10 * t3 - t12 * t2;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t13 * t6, t6, t8 * t7, t7 * t2, 0, 0; t11 * t6 + t13 * t7, -t7, t8 * t6, t6 * t2, 0, 0; 0, 0, t9, -t3, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (87->25), mult. (98->35), div. (0->0), fcn. (111->7), ass. (0->22)
	t18 = pkin(3) + pkin(8) + r_i_i_C(3);
	t8 = pkin(9) + qJ(3);
	t7 = cos(t8);
	t16 = t18 * t7;
	t6 = sin(t8);
	t24 = t16 + t6 * qJ(4) + cos(pkin(9)) * pkin(2) + pkin(1);
	t23 = pkin(4) + pkin(7) + qJ(2);
	t10 = sin(qJ(5));
	t13 = cos(qJ(1));
	t22 = t10 * t13;
	t11 = sin(qJ(1));
	t21 = t11 * t10;
	t12 = cos(qJ(5));
	t20 = t11 * t12;
	t19 = t12 * t13;
	t15 = r_i_i_C(1) * t10 + r_i_i_C(2) * t12 + qJ(4);
	t14 = t15 * t7 - t18 * t6;
	t4 = -t6 * t21 + t19;
	t3 = t6 * t20 + t22;
	t2 = t6 * t22 + t20;
	t1 = t6 * t19 - t21;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t24 * t11 + t23 * t13, t11, t14 * t13, t13 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t23 * t11 + t24 * t13, -t13, t14 * t11, t11 * t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 0, t15 * t6 + t16, -t7, (-r_i_i_C(1) * t12 + r_i_i_C(2) * t10) * t7, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:37:41
	% EndTime: 2019-10-10 00:37:41
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (115->29), mult. (128->39), div. (0->0), fcn. (144->7), ass. (0->23)
	t12 = sin(qJ(5));
	t20 = pkin(3) + r_i_i_C(3) + qJ(6) + pkin(8);
	t9 = pkin(9) + qJ(3);
	t8 = cos(t9);
	t18 = t20 * t8;
	t7 = sin(t9);
	t28 = -(pkin(5) * t12 + qJ(4)) * t7 - cos(pkin(9)) * pkin(2) - pkin(1) - t18;
	t26 = pkin(5) + r_i_i_C(1);
	t14 = cos(qJ(5));
	t25 = pkin(5) * t14 + pkin(4) + pkin(7) + qJ(2);
	t15 = cos(qJ(1));
	t24 = t12 * t15;
	t13 = sin(qJ(1));
	t23 = t13 * t12;
	t22 = t13 * t14;
	t21 = t14 * t15;
	t1 = t7 * t21 - t23;
	t3 = t7 * t22 + t24;
	t17 = r_i_i_C(2) * t14 + t26 * t12 + qJ(4);
	t16 = t17 * t8 - t20 * t7;
	t4 = -t7 * t23 + t21;
	t2 = t7 * t24 + t22;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t28 * t13 + t25 * t15, t13, t16 * t15, t15 * t7, -t2 * r_i_i_C(2) + t26 * t1, t15 * t8; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t25 * t13 - t28 * t15, -t15, t16 * t13, t13 * t7, t4 * r_i_i_C(2) + t26 * t3, t13 * t8; 0, 0, t17 * t7 + t18, -t8, (r_i_i_C(2) * t12 - t26 * t14) * t8, t7;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end