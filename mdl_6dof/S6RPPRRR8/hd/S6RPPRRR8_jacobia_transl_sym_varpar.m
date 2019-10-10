% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRR8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
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
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (11->6), mult. (16->6), div. (0->0), fcn. (20->4), ass. (0->5)
	t6 = pkin(1) + r_i_i_C(3) + qJ(3);
	t5 = r_i_i_C(1) * sin(pkin(10)) + r_i_i_C(2) * cos(pkin(10)) + qJ(2);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t6 * t3 + t5 * t4, t3, t4, 0, 0, 0; t5 * t3 + t6 * t4, -t4, t3, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (28->10), mult. (30->11), div. (0->0), fcn. (34->5), ass. (0->10)
	t11 = pkin(1) + r_i_i_C(3) + pkin(7) + qJ(3);
	t3 = pkin(10) + qJ(4);
	t1 = sin(t3);
	t2 = cos(t3);
	t10 = r_i_i_C(1) * t2 - r_i_i_C(2) * t1;
	t9 = -t1 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t8 = pkin(3) * sin(pkin(10)) + qJ(2) - t9;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t4 = [-t11 * t6 + t8 * t7, t6, t7, t10 * t6, 0, 0; t11 * t7 + t8 * t6, -t7, t6, -t10 * t7, 0, 0; 0, 0, 0, t9, 0, 0;];
	Ja_transl = t4;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (73->25), mult. (89->33), div. (0->0), fcn. (101->7), ass. (0->22)
	t20 = pkin(8) + r_i_i_C(3);
	t7 = pkin(10) + qJ(4);
	t6 = cos(t7);
	t23 = t20 * t6;
	t10 = sin(qJ(5));
	t12 = cos(qJ(5));
	t15 = r_i_i_C(1) * t12 - r_i_i_C(2) * t10 + pkin(4);
	t5 = sin(t7);
	t22 = t15 * t6 + t20 * t5;
	t21 = pkin(1) + pkin(7) + qJ(3);
	t13 = cos(qJ(1));
	t19 = t10 * t13;
	t11 = sin(qJ(1));
	t18 = t11 * t10;
	t17 = t11 * t12;
	t16 = t12 * t13;
	t14 = pkin(3) * sin(pkin(10)) + t5 * pkin(4) - t23 + qJ(2);
	t4 = t5 * t16 - t18;
	t3 = t5 * t19 + t17;
	t2 = t5 * t17 + t19;
	t1 = -t5 * t18 + t16;
	t8 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t21 * t11 + t14 * t13, t11, t13, t22 * t11, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t14 * t11 + t21 * t13, -t13, t11, -t22 * t13, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 0, 0, -t15 * t5 + t23, (-r_i_i_C(1) * t10 - r_i_i_C(2) * t12) * t6, 0;];
	Ja_transl = t8;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:13:28
	% EndTime: 2019-10-10 00:13:28
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (140->32), mult. (130->43), div. (0->0), fcn. (146->9), ass. (0->31)
	t14 = pkin(10) + qJ(4);
	t11 = cos(t14);
	t32 = r_i_i_C(3) + pkin(9) + pkin(8);
	t37 = t32 * t11;
	t10 = sin(t14);
	t15 = qJ(5) + qJ(6);
	t12 = sin(t15);
	t13 = cos(t15);
	t20 = cos(qJ(5));
	t9 = pkin(5) * t20 + pkin(4);
	t24 = r_i_i_C(1) * t13 - r_i_i_C(2) * t12 + t9;
	t36 = t32 * t10 + t24 * t11;
	t21 = cos(qJ(1));
	t27 = t13 * t21;
	t19 = sin(qJ(1));
	t30 = t12 * t19;
	t5 = -t10 * t30 + t27;
	t28 = t13 * t19;
	t29 = t12 * t21;
	t6 = t10 * t28 + t29;
	t35 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t7 = t10 * t29 + t28;
	t8 = t10 * t27 - t30;
	t34 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t18 = sin(qJ(5));
	t33 = pkin(5) * t18;
	t31 = t10 * t18;
	t26 = pkin(1) + pkin(7) + qJ(3) + t33;
	t25 = -r_i_i_C(1) * t12 - r_i_i_C(2) * t13;
	t23 = pkin(3) * sin(pkin(10)) + t10 * t9 - t37 + qJ(2);
	t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t26 * t19 + t23 * t21, t19, t21, t36 * t19, (-t19 * t31 + t20 * t21) * pkin(5) + t35, t35; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t23 * t19 + t26 * t21, -t21, t19, -t36 * t21, (t19 * t20 + t21 * t31) * pkin(5) + t34, t34; 0, 0, 0, -t24 * t10 + t37, (t25 - t33) * t11, t25 * t11;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end