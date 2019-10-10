% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:37
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:30
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:30
	% EndTime: 2019-10-10 09:37:31
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
	% StartTime: 2019-10-10 09:37:30
	% EndTime: 2019-10-10 09:37:31
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
	% StartTime: 2019-10-10 09:37:30
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (29->10), mult. (31->12), div. (0->0), fcn. (33->6), ass. (0->10)
	t5 = qJ(2) + pkin(10);
	t2 = sin(t5);
	t3 = cos(t5);
	t14 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2) + cos(qJ(2)) * pkin(2);
	t13 = r_i_i_C(3) + qJ(3) + pkin(7);
	t11 = pkin(1) + t14;
	t10 = -sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t11 * t8 + t13 * t9, t10 * t9, t8, 0, 0, 0; t11 * t9 + t13 * t8, t10 * t8, -t9, 0, 0, 0; 0, t14, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:30
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (65->16), mult. (72->18), div. (0->0), fcn. (81->8), ass. (0->14)
	t6 = sin(pkin(11));
	t7 = cos(pkin(11));
	t15 = r_i_i_C(1) * t7 - r_i_i_C(2) * t6 + pkin(3);
	t16 = r_i_i_C(3) + qJ(4);
	t5 = qJ(2) + pkin(10);
	t2 = sin(t5);
	t3 = cos(t5);
	t18 = cos(qJ(2)) * pkin(2) + t15 * t3 + t16 * t2;
	t17 = pkin(1) + t18;
	t14 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7 + pkin(7) + qJ(3);
	t12 = -sin(qJ(2)) * pkin(2) - t15 * t2 + t16 * t3;
	t11 = cos(qJ(1));
	t10 = sin(qJ(1));
	t1 = [-t17 * t10 + t14 * t11, t12 * t11, t10, t11 * t2, 0, 0; t14 * t10 + t17 * t11, t12 * t10, -t11, t10 * t2, 0, 0; 0, t18, 0, -t3, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:30
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (114->29), mult. (101->40), div. (0->0), fcn. (114->10), ass. (0->22)
	t27 = r_i_i_C(3) + pkin(8) + qJ(4);
	t13 = qJ(2) + pkin(10);
	t8 = sin(t13);
	t29 = cos(qJ(2)) * pkin(2) + t27 * t8;
	t10 = cos(t13);
	t5 = cos(pkin(11)) * pkin(4) + pkin(3);
	t28 = t10 * t5 + pkin(1) + t29;
	t18 = sin(qJ(1));
	t26 = t10 * t18;
	t19 = cos(qJ(1));
	t25 = t10 * t19;
	t22 = sin(pkin(11)) * pkin(4) + qJ(3) + pkin(7);
	t12 = pkin(11) + qJ(5);
	t7 = sin(t12);
	t9 = cos(t12);
	t21 = r_i_i_C(1) * t9 - r_i_i_C(2) * t7 + t5;
	t20 = -sin(qJ(2)) * pkin(2) + t27 * t10 - t21 * t8;
	t4 = t18 * t7 + t9 * t25;
	t3 = t18 * t9 - t7 * t25;
	t2 = t19 * t7 - t9 * t26;
	t1 = t19 * t9 + t7 * t26;
	t6 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t28 * t18 + t22 * t19, t20 * t19, t18, t19 * t8, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t22 * t18 + t28 * t19, t20 * t18, -t19, t18 * t8, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, t21 * t10 + t29, 0, -t10, (-r_i_i_C(1) * t7 - r_i_i_C(2) * t9) * t8, 0;];
	Ja_transl = t6;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:37:30
	% EndTime: 2019-10-10 09:37:31
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (203->37), mult. (140->49), div. (0->0), fcn. (157->12), ass. (0->30)
	t22 = qJ(2) + pkin(10);
	t15 = sin(t22);
	t36 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
	t41 = cos(qJ(2)) * pkin(2) + t36 * t15;
	t17 = cos(t22);
	t21 = pkin(11) + qJ(5);
	t16 = cos(t21);
	t9 = pkin(5) * t16 + cos(pkin(11)) * pkin(4) + pkin(3);
	t40 = t17 * t9 + pkin(1) + t41;
	t18 = qJ(6) + t21;
	t12 = cos(t18);
	t26 = cos(qJ(1));
	t11 = sin(t18);
	t25 = sin(qJ(1));
	t34 = t25 * t11;
	t5 = t12 * t26 + t17 * t34;
	t33 = t25 * t12;
	t6 = t11 * t26 - t17 * t33;
	t39 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t35 = t17 * t26;
	t7 = -t11 * t35 + t33;
	t8 = t12 * t35 + t34;
	t38 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t14 = sin(t21);
	t37 = pkin(5) * t14;
	t32 = t37 + sin(pkin(11)) * pkin(4) + qJ(3) + pkin(7);
	t29 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
	t28 = r_i_i_C(1) * t12 - r_i_i_C(2) * t11 + t9;
	t27 = -sin(qJ(2)) * pkin(2) - t28 * t15 + t36 * t17;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t40 * t25 + t32 * t26, t27 * t26, t25, t26 * t15, (-t14 * t35 + t16 * t25) * pkin(5) + t38, t38; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t32 * t25 + t40 * t26, t27 * t25, -t26, t25 * t15, (-t14 * t17 * t25 - t16 * t26) * pkin(5) + t39, t39; 0, t17 * t28 + t41, 0, -t17, (t29 - t37) * t15, t29 * t15;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end