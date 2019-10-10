% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:02
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR10_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR10_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR10_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
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
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
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
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (14->8), mult. (26->10), div. (0->0), fcn. (28->4), ass. (0->9)
	t8 = pkin(1) + pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(3));
	t3 = cos(qJ(3));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t5 = qJ(2) - t6;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t8 * t2 + t5 * t4, t2, t7 * t2, 0, 0, 0; t5 * t2 + t8 * t4, -t4, -t7 * t4, 0, 0, 0; 0, 0, t6, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (32->14), mult. (67->16), div. (0->0), fcn. (76->6), ass. (0->13)
	t10 = r_i_i_C(3) + qJ(4);
	t3 = sin(qJ(3));
	t5 = cos(qJ(3));
	t1 = sin(pkin(10));
	t2 = cos(pkin(10));
	t9 = r_i_i_C(1) * t2 - r_i_i_C(2) * t1 + pkin(3);
	t12 = t10 * t5 - t9 * t3;
	t11 = t10 * t3 + t9 * t5;
	t8 = r_i_i_C(1) * t1 + r_i_i_C(2) * t2 + pkin(1) + pkin(7);
	t7 = qJ(2) - t12;
	t6 = cos(qJ(1));
	t4 = sin(qJ(1));
	t13 = [-t8 * t4 + t7 * t6, t4, t11 * t4, -t4 * t5, 0, 0; t7 * t4 + t8 * t6, -t6, -t11 * t6, t6 * t5, 0, 0; 0, 0, t12, t3, 0, 0;];
	Ja_transl = t13;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (75->27), mult. (96->38), div. (0->0), fcn. (109->8), ass. (0->21)
	t13 = cos(qJ(3));
	t20 = r_i_i_C(3) + pkin(8) + qJ(4);
	t22 = t20 * t13;
	t11 = sin(qJ(3));
	t5 = cos(pkin(10)) * pkin(4) + pkin(3);
	t8 = pkin(10) + qJ(5);
	t6 = sin(t8);
	t7 = cos(t8);
	t16 = r_i_i_C(1) * t7 - r_i_i_C(2) * t6 + t5;
	t21 = t20 * t11 + t16 * t13;
	t12 = sin(qJ(1));
	t19 = t11 * t12;
	t14 = cos(qJ(1));
	t18 = t11 * t14;
	t17 = sin(pkin(10)) * pkin(4) + pkin(1) + pkin(7);
	t15 = t11 * t5 + qJ(2) - t22;
	t4 = -t12 * t6 + t7 * t18;
	t3 = t12 * t7 + t6 * t18;
	t2 = t14 * t6 + t7 * t19;
	t1 = t14 * t7 - t6 * t19;
	t9 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t17 * t12 + t15 * t14, t12, t21 * t12, -t12 * t13, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * t12 + t17 * t14, -t14, -t21 * t14, t14 * t13, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 0, -t16 * t11 + t22, t11, (-r_i_i_C(1) * t6 - r_i_i_C(2) * t7) * t13, 0;];
	Ja_transl = t9;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:02:11
	% EndTime: 2019-10-10 01:02:11
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (155->35), mult. (135->47), div. (0->0), fcn. (152->10), ass. (0->29)
	t20 = cos(qJ(3));
	t29 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
	t34 = t29 * t20;
	t18 = sin(qJ(3));
	t17 = pkin(10) + qJ(5);
	t15 = qJ(6) + t17;
	t11 = sin(t15);
	t12 = cos(t15);
	t14 = cos(t17);
	t9 = pkin(5) * t14 + cos(pkin(10)) * pkin(4) + pkin(3);
	t23 = r_i_i_C(1) * t12 - r_i_i_C(2) * t11 + t9;
	t33 = t29 * t18 + t23 * t20;
	t21 = cos(qJ(1));
	t19 = sin(qJ(1));
	t27 = t19 * t11;
	t5 = t12 * t21 - t18 * t27;
	t26 = t19 * t12;
	t6 = t11 * t21 + t18 * t26;
	t32 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t28 = t18 * t21;
	t7 = t11 * t28 + t26;
	t8 = t12 * t28 - t27;
	t31 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t13 = sin(t17);
	t30 = pkin(5) * t13;
	t25 = pkin(1) + pkin(7) + t30 + sin(pkin(10)) * pkin(4);
	t24 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
	t22 = t18 * t9 + qJ(2) - t34;
	t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t25 * t19 + t22 * t21, t19, t33 * t19, -t19 * t20, (-t13 * t18 * t19 + t14 * t21) * pkin(5) + t32, t32; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t22 * t19 + t25 * t21, -t21, -t33 * t21, t21 * t20, (t13 * t28 + t14 * t19) * pkin(5) + t31, t31; 0, 0, -t23 * t18 + t34, t18, (t24 - t30) * t20, t24 * t20;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end