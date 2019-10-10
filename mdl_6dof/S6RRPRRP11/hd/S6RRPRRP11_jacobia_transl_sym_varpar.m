% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP11_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP11_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP11_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
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
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t5 * t2 + t8 * t4, t6 * t4, 0, 0, 0, 0; t8 * t2 + t5 * t4, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
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
	% StartTime: 2019-10-10 10:46:42
	% EndTime: 2019-10-10 10:46:43
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
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (101->30), mult. (144->43), div. (0->0), fcn. (159->8), ass. (0->28)
	t15 = sin(qJ(2));
	t18 = cos(qJ(2));
	t25 = pkin(2) + r_i_i_C(3) + pkin(9) + pkin(8);
	t23 = t25 * t18;
	t14 = sin(qJ(4));
	t24 = pkin(4) * t14 + qJ(3);
	t34 = -t24 * t15 - pkin(1) - t23;
	t13 = qJ(4) + qJ(5);
	t11 = sin(t13);
	t12 = cos(t13);
	t16 = sin(qJ(1));
	t19 = cos(qJ(1));
	t26 = t19 * t15;
	t5 = -t11 * t16 + t12 * t26;
	t6 = t11 * t26 + t12 * t16;
	t32 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t27 = t16 * t15;
	t7 = t11 * t19 + t12 * t27;
	t8 = -t11 * t27 + t12 * t19;
	t31 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t17 = cos(qJ(4));
	t30 = pkin(4) * t17;
	t29 = r_i_i_C(1) * t12;
	t28 = pkin(7) + pkin(3) + t30;
	t22 = r_i_i_C(1) * t11 + r_i_i_C(2) * t12 + t24;
	t21 = -t25 * t15 + t22 * t18;
	t9 = t18 * t11 * r_i_i_C(2);
	t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t34 * t16 + t28 * t19, t21 * t19, t26, (-t14 * t16 + t17 * t26) * pkin(4) + t32, t32, 0; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t28 * t16 - t34 * t19, t21 * t16, t27, (t14 * t19 + t17 * t27) * pkin(4) + t31, t31, 0; 0, t22 * t15 + t23, -t18, t9 + (-t29 - t30) * t18, -t18 * t29 + t9, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:46:43
	% EndTime: 2019-10-10 10:46:43
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (140->35), mult. (170->48), div. (0->0), fcn. (188->8), ass. (0->26)
	t19 = sin(qJ(2));
	t21 = cos(qJ(2));
	t26 = pkin(2) + r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
	t25 = t26 * t21;
	t18 = qJ(4) + qJ(5);
	t14 = sin(t18);
	t10 = pkin(5) * t14 + sin(qJ(4)) * pkin(4);
	t27 = qJ(3) + t10;
	t34 = -t27 * t19 - pkin(1) - t25;
	t15 = cos(t18);
	t11 = pkin(5) * t15 + cos(qJ(4)) * pkin(4);
	t32 = pkin(7) + pkin(3) + t11;
	t20 = sin(qJ(1));
	t22 = cos(qJ(1));
	t28 = t22 * t19;
	t5 = -t20 * t14 + t15 * t28;
	t6 = t14 * t28 + t20 * t15;
	t31 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t29 = t20 * t19;
	t7 = t22 * t14 + t15 * t29;
	t8 = -t14 * t29 + t22 * t15;
	t30 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t24 = r_i_i_C(1) * t14 + r_i_i_C(2) * t15 + t27;
	t23 = -t26 * t19 + t24 * t21;
	t12 = t21 * t14 * r_i_i_C(2);
	t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t34 * t20 + t32 * t22, t23 * t22, t28, -t20 * t10 + t11 * t28 + t31, t5 * pkin(5) + t31, t22 * t21; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t32 * t20 - t34 * t22, t23 * t20, t29, t22 * t10 + t11 * t29 + t30, t7 * pkin(5) + t30, t20 * t21; 0, t24 * t19 + t25, -t21, t12 + (-r_i_i_C(1) * t15 - t11) * t21, t12 + (-pkin(5) - r_i_i_C(1)) * t21 * t15, t19;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end