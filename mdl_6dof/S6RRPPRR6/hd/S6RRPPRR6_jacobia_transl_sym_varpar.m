% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
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
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
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
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (22->9), mult. (44->12), div. (0->0), fcn. (47->4), ass. (0->11)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(3) + qJ(3);
	t9 = pkin(2) + r_i_i_C(1);
	t6 = t7 * t1 + t3 * t9;
	t10 = pkin(1) + t6;
	t8 = pkin(7) + r_i_i_C(2);
	t5 = -t1 * t9 + t7 * t3;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t11 = [-t10 * t2 + t8 * t4, t5 * t4, t4 * t1, 0, 0, 0; t10 * t4 + t8 * t2, t5 * t2, t2 * t1, 0, 0, 0; 0, t6, -t3, 0, 0, 0;];
	Ja_transl = t11;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (40->24), mult. (83->34), div. (0->0), fcn. (98->6), ass. (0->17)
	t16 = pkin(2) + pkin(3);
	t7 = sin(qJ(2));
	t9 = cos(qJ(2));
	t17 = t7 * qJ(3) + t16 * t9 + pkin(1);
	t15 = pkin(7) - r_i_i_C(3) - qJ(4);
	t5 = sin(pkin(10));
	t6 = cos(pkin(10));
	t13 = t5 * t9 - t6 * t7;
	t12 = t5 * t7 + t6 * t9;
	t11 = qJ(3) * t9 - t16 * t7;
	t10 = cos(qJ(1));
	t8 = sin(qJ(1));
	t4 = t12 * t10;
	t3 = t13 * t10;
	t2 = t12 * t8;
	t1 = t13 * t8;
	t14 = [-t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * t10 - t17 * t8, t3 * r_i_i_C(1) + t4 * r_i_i_C(2) + t11 * t10, t10 * t7, -t8, 0, 0; t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t17 * t10 + t15 * t8, t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + t11 * t8, t8 * t7, t10, 0, 0; 0, (r_i_i_C(1) * t5 + r_i_i_C(2) * t6 + qJ(3)) * t7 + (r_i_i_C(1) * t6 - r_i_i_C(2) * t5 + t16) * t9, -t9, 0, 0, 0;];
	Ja_transl = t14;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (93->25), mult. (128->33), div. (0->0), fcn. (149->8), ass. (0->25)
	t14 = sin(qJ(2));
	t16 = cos(qJ(2));
	t24 = pkin(4) * sin(pkin(10)) + qJ(3);
	t28 = pkin(2) + cos(pkin(10)) * pkin(4) + pkin(3);
	t19 = t24 * t14 + t28 * t16;
	t29 = pkin(1) + t19;
	t11 = pkin(10) + qJ(5);
	t9 = sin(t11);
	t27 = t16 * t9;
	t17 = cos(qJ(1));
	t26 = t17 * t14;
	t25 = pkin(7) - r_i_i_C(3) - pkin(8) - qJ(4);
	t15 = sin(qJ(1));
	t10 = cos(t11);
	t6 = t10 * t14 - t27;
	t1 = t6 * t15;
	t20 = t10 * t16 + t14 * t9;
	t2 = t20 * t15;
	t23 = r_i_i_C(1) * t1 - r_i_i_C(2) * t2;
	t3 = -t10 * t26 + t17 * t27;
	t4 = t20 * t17;
	t22 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t21 = -r_i_i_C(1) * t20 - t6 * r_i_i_C(2);
	t18 = -t28 * t14 + t24 * t16;
	t5 = [-t2 * r_i_i_C(1) - t1 * r_i_i_C(2) - t29 * t15 + t25 * t17, t18 * t17 - t22, t26, -t15, t22, 0; t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t25 * t15 + t29 * t17, t18 * t15 - t23, t15 * t14, t17, t23, 0; 0, t19 - t21, -t16, 0, t21, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:44:54
	% EndTime: 2019-10-10 09:44:54
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (216->39), mult. (290->55), div. (0->0), fcn. (349->10), ass. (0->32)
	t18 = sin(qJ(2));
	t21 = cos(qJ(2));
	t33 = pkin(4) * sin(pkin(10)) + qJ(3);
	t36 = pkin(2) + cos(pkin(10)) * pkin(4) + pkin(3);
	t24 = t33 * t18 + t36 * t21;
	t41 = pkin(1) + t24;
	t34 = pkin(10) + qJ(5);
	t31 = sin(t34);
	t32 = cos(t34);
	t8 = t18 * t32 - t21 * t31;
	t17 = sin(qJ(6));
	t20 = cos(qJ(6));
	t25 = t20 * r_i_i_C(1) - t17 * r_i_i_C(2) + pkin(5);
	t19 = sin(qJ(1));
	t3 = t8 * t19;
	t37 = -r_i_i_C(3) - pkin(9);
	t7 = t18 * t31 + t21 * t32;
	t4 = t7 * t19;
	t40 = t25 * t3 - t37 * t4;
	t22 = cos(qJ(1));
	t26 = t22 * t31;
	t27 = t22 * t32;
	t5 = -t18 * t26 - t21 * t27;
	t6 = -t18 * t27 + t21 * t26;
	t39 = -t25 * t6 + t37 * t5;
	t38 = -t25 * t7 - t37 * t8;
	t35 = pkin(7) - pkin(8) - qJ(4);
	t30 = -t17 * r_i_i_C(1) - t20 * r_i_i_C(2);
	t23 = -t36 * t18 + t33 * t21;
	t2 = -t19 * t17 - t5 * t20;
	t1 = t5 * t17 - t19 * t20;
	t9 = [-t37 * t3 - t25 * t4 + (t30 + t35) * t22 - t41 * t19, t23 * t22 - t39, t22 * t18, -t19, t39, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -t5 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t35 * t19 + t41 * t22 - t37 * t6, t23 * t19 - t40, t19 * t18, t22, t40, (-t4 * t17 + t22 * t20) * r_i_i_C(1) + (-t22 * t17 - t4 * t20) * r_i_i_C(2); 0, t24 - t38, -t21, 0, t38, t30 * t8;];
	Ja_transl = t9;
else
	Ja_transl=NaN(3,6);
end