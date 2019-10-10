% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP9_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP9_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP9_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
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
	% StartTime: 2019-10-10 08:52:21
	% EndTime: 2019-10-10 08:52:22
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
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (38->23), mult. (85->32), div. (0->0), fcn. (95->6), ass. (0->21)
	t15 = pkin(8) + r_i_i_C(3);
	t9 = cos(qJ(3));
	t20 = t15 * t9;
	t5 = sin(qJ(4));
	t8 = cos(qJ(4));
	t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5 + pkin(3);
	t6 = sin(qJ(3));
	t19 = t12 * t9 + t15 * t6;
	t7 = sin(qJ(1));
	t18 = t7 * t5;
	t17 = t7 * t8;
	t16 = pkin(1) + pkin(7);
	t10 = cos(qJ(1));
	t14 = t10 * t5;
	t13 = t10 * t8;
	t11 = t6 * pkin(3) + qJ(2) - t20;
	t4 = t6 * t13 - t18;
	t3 = t6 * t14 + t17;
	t2 = t6 * t17 + t14;
	t1 = -t6 * t18 + t13;
	t21 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t11 * t10 - t16 * t7, t7, t19 * t7, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t10 + t11 * t7, -t10, -t19 * t10, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0, 0; 0, 0, -t12 * t6 + t20, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t9, 0, 0;];
	Ja_transl = t21;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (96->30), mult. (126->43), div. (0->0), fcn. (140->8), ass. (0->27)
	t17 = cos(qJ(3));
	t26 = r_i_i_C(3) + pkin(9) + pkin(8);
	t31 = t26 * t17;
	t14 = sin(qJ(3));
	t12 = qJ(4) + qJ(5);
	t10 = sin(t12);
	t11 = cos(t12);
	t16 = cos(qJ(4));
	t9 = pkin(4) * t16 + pkin(3);
	t21 = r_i_i_C(1) * t11 - r_i_i_C(2) * t10 + t9;
	t30 = t26 * t14 + t21 * t17;
	t18 = cos(qJ(1));
	t15 = sin(qJ(1));
	t25 = t14 * t15;
	t5 = -t10 * t25 + t11 * t18;
	t6 = t10 * t18 + t11 * t25;
	t29 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t24 = t14 * t18;
	t7 = t10 * t24 + t11 * t15;
	t8 = -t10 * t15 + t11 * t24;
	t28 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t13 = sin(qJ(4));
	t27 = pkin(4) * t13;
	t23 = pkin(1) + pkin(7) + t27;
	t22 = -r_i_i_C(1) * t10 - r_i_i_C(2) * t11;
	t20 = t14 * t9 + qJ(2) - t31;
	t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t23 * t15 + t20 * t18, t15, t30 * t15, (-t13 * t25 + t16 * t18) * pkin(4) + t29, t29, 0; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t20 * t15 + t23 * t18, -t18, -t30 * t18, (t13 * t24 + t15 * t16) * pkin(4) + t28, t28, 0; 0, 0, -t21 * t14 + t31, (t22 - t27) * t17, t22 * t17, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:52:22
	% EndTime: 2019-10-10 08:52:22
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (136->37), mult. (152->48), div. (0->0), fcn. (169->8), ass. (0->27)
	t20 = cos(qJ(3));
	t28 = r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
	t33 = t28 * t20;
	t18 = sin(qJ(3));
	t17 = qJ(4) + qJ(5);
	t13 = sin(t17);
	t14 = cos(t17);
	t11 = pkin(5) * t14 + cos(qJ(4)) * pkin(4);
	t9 = pkin(3) + t11;
	t23 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
	t32 = t28 * t18 + t23 * t20;
	t21 = cos(qJ(1));
	t19 = sin(qJ(1));
	t26 = t19 * t13;
	t5 = t14 * t21 - t18 * t26;
	t25 = t19 * t14;
	t6 = t13 * t21 + t18 * t25;
	t31 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t27 = t18 * t21;
	t7 = t13 * t27 + t25;
	t8 = t14 * t27 - t26;
	t30 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t29 = r_i_i_C(2) * t14;
	t10 = pkin(5) * t13 + sin(qJ(4)) * pkin(4);
	t24 = pkin(1) + pkin(7) + t10;
	t22 = t18 * t9 + qJ(2) - t33;
	t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) - t24 * t19 + t22 * t21, t19, t32 * t19, -t19 * t18 * t10 + t11 * t21 + t31, t5 * pkin(5) + t31, -t19 * t20; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t22 * t19 + t24 * t21, -t21, -t32 * t21, t10 * t27 + t19 * t11 + t30, t7 * pkin(5) + t30, t21 * t20; 0, 0, -t23 * t18 + t33, (-r_i_i_C(1) * t13 - t10 - t29) * t20, (-t29 + (-pkin(5) - r_i_i_C(1)) * t13) * t20, t18;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end