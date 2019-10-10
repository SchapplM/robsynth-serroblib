% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:22
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR11_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR11_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR11_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:42
	% EndTime: 2019-10-10 10:22:42
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:42
	% EndTime: 2019-10-10 10:22:42
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
	% StartTime: 2019-10-10 10:22:42
	% EndTime: 2019-10-10 10:22:42
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
	% StartTime: 2019-10-10 10:22:42
	% EndTime: 2019-10-10 10:22:42
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
	% StartTime: 2019-10-10 10:22:42
	% EndTime: 2019-10-10 10:22:43
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (43->22), mult. (96->34), div. (0->0), fcn. (107->6), ass. (0->21)
	t15 = pkin(2) + pkin(8) + r_i_i_C(3);
	t9 = cos(qJ(2));
	t13 = t15 * t9;
	t6 = sin(qJ(2));
	t21 = t6 * qJ(3) + pkin(1) + t13;
	t7 = sin(qJ(1));
	t20 = t7 * t6;
	t8 = cos(qJ(4));
	t19 = t7 * t8;
	t18 = pkin(3) + pkin(7);
	t10 = cos(qJ(1));
	t17 = t10 * t6;
	t16 = t10 * t8;
	t5 = sin(qJ(4));
	t12 = r_i_i_C(1) * t5 + r_i_i_C(2) * t8 + qJ(3);
	t11 = t12 * t9 - t15 * t6;
	t4 = -t5 * t20 + t16;
	t3 = t10 * t5 + t6 * t19;
	t2 = t5 * t17 + t19;
	t1 = t6 * t16 - t7 * t5;
	t14 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t18 * t10 - t21 * t7, t11 * t10, t17, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t21 * t10 + t18 * t7, t11 * t7, t20, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0, 0; 0, t12 * t6 + t13, -t9, (-r_i_i_C(1) * t8 + r_i_i_C(2) * t5) * t9, 0, 0;];
	Ja_transl = t14;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:42
	% EndTime: 2019-10-10 10:22:43
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (84->29), mult. (126->44), div. (0->0), fcn. (140->8), ass. (0->24)
	t11 = sin(qJ(2));
	t14 = cos(qJ(2));
	t20 = pkin(2) + r_i_i_C(3) + qJ(5) + pkin(8);
	t18 = t20 * t14;
	t10 = sin(qJ(4));
	t19 = pkin(4) * t10 + qJ(3);
	t26 = -t19 * t11 - pkin(1) - t18;
	t13 = cos(qJ(4));
	t23 = pkin(4) * t13;
	t24 = pkin(7) + pkin(3) + t23;
	t12 = sin(qJ(1));
	t22 = t12 * t11;
	t15 = cos(qJ(1));
	t21 = t15 * t11;
	t8 = qJ(4) + pkin(10);
	t6 = sin(t8);
	t7 = cos(t8);
	t17 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7 + t19;
	t16 = -t20 * t11 + t17 * t14;
	t4 = t15 * t7 - t6 * t22;
	t3 = t15 * t6 + t7 * t22;
	t2 = t12 * t7 + t6 * t21;
	t1 = -t12 * t6 + t7 * t21;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t26 * t12 + t24 * t15, t16 * t15, t21, t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + (-t10 * t12 + t13 * t21) * pkin(4), t15 * t14, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t24 * t12 - t26 * t15, t16 * t12, t22, t3 * r_i_i_C(1) + t4 * r_i_i_C(2) + (t10 * t15 + t13 * t22) * pkin(4), t12 * t14, 0; 0, t17 * t11 + t18, -t14, (-r_i_i_C(1) * t7 + r_i_i_C(2) * t6 - t23) * t14, t11, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:22:42
	% EndTime: 2019-10-10 10:22:43
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (164->33), mult. (158->45), div. (0->0), fcn. (176->10), ass. (0->28)
	t20 = sin(qJ(2));
	t22 = cos(qJ(2));
	t27 = pkin(2) + r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
	t26 = t27 * t22;
	t19 = qJ(4) + pkin(10);
	t10 = pkin(5) * sin(t19) + sin(qJ(4)) * pkin(4);
	t28 = qJ(3) + t10;
	t36 = -t28 * t20 - pkin(1) - t26;
	t11 = pkin(5) * cos(t19) + cos(qJ(4)) * pkin(4);
	t34 = pkin(7) + pkin(3) + t11;
	t16 = qJ(6) + t19;
	t14 = sin(t16);
	t15 = cos(t16);
	t21 = sin(qJ(1));
	t23 = cos(qJ(1));
	t29 = t23 * t20;
	t5 = -t21 * t14 + t15 * t29;
	t6 = t14 * t29 + t21 * t15;
	t33 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t30 = t21 * t20;
	t7 = t23 * t14 + t15 * t30;
	t8 = -t14 * t30 + t23 * t15;
	t32 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t31 = r_i_i_C(1) * t15;
	t25 = r_i_i_C(1) * t14 + r_i_i_C(2) * t15 + t28;
	t24 = -t27 * t20 + t25 * t22;
	t12 = t22 * t14 * r_i_i_C(2);
	t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t36 * t21 + t34 * t23, t24 * t23, t29, -t21 * t10 + t11 * t29 + t33, t23 * t22, t33; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t34 * t21 - t36 * t23, t24 * t21, t30, t23 * t10 + t11 * t30 + t32, t21 * t22, t32; 0, t25 * t20 + t26, -t22, t12 + (-t11 - t31) * t22, t20, -t22 * t31 + t12;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end