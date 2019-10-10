% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
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
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(10);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (27->9), mult. (26->12), div. (0->0), fcn. (26->6), ass. (0->10)
	t9 = pkin(7) + r_i_i_C(3);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t7 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t6 = pkin(2) + t8;
	t3 = qJ(1) + pkin(10);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [-sin(qJ(1)) * pkin(1) + t9 * t2 - t6 * t1, 0, t7 * t2, 0, 0, 0; cos(qJ(1)) * pkin(1) + t9 * t1 + t6 * t2, 0, t7 * t1, 0, 0, 0; 0, 1, t8, 0, 0, 0;];
	Ja_transl = t10;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (73->23), mult. (85->36), div. (0->0), fcn. (93->8), ass. (0->19)
	t11 = cos(qJ(3));
	t18 = pkin(8) + r_i_i_C(3);
	t9 = sin(qJ(3));
	t15 = t18 * t9;
	t19 = t11 * pkin(3) + pkin(2) + t15;
	t8 = sin(qJ(4));
	t17 = t11 * t8;
	t10 = cos(qJ(4));
	t16 = t10 * t11;
	t13 = r_i_i_C(1) * t10 - r_i_i_C(2) * t8 + pkin(3);
	t12 = t18 * t11 - t13 * t9;
	t7 = qJ(1) + pkin(10);
	t6 = cos(t7);
	t5 = sin(t7);
	t4 = t6 * t16 + t5 * t8;
	t3 = t5 * t10 - t6 * t17;
	t2 = -t5 * t16 + t6 * t8;
	t1 = t6 * t10 + t5 * t17;
	t14 = [-sin(qJ(1)) * pkin(1) + t6 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t19 * t5, 0, t12 * t6, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; cos(qJ(1)) * pkin(1) + t5 * pkin(7) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t19 * t6, 0, t12 * t5, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, 1, t13 * t11 + t15, (-r_i_i_C(1) * t8 - r_i_i_C(2) * t10) * t9, 0, 0;];
	Ja_transl = t14;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (145->31), mult. (126->46), div. (0->0), fcn. (138->10), ass. (0->29)
	t19 = cos(qJ(3));
	t17 = sin(qJ(3));
	t30 = r_i_i_C(3) + pkin(9) + pkin(8);
	t24 = t30 * t17;
	t18 = cos(qJ(4));
	t9 = pkin(4) * t18 + pkin(3);
	t34 = t19 * t9 + pkin(2) + t24;
	t14 = qJ(1) + pkin(10);
	t10 = sin(t14);
	t11 = cos(t14);
	t15 = qJ(4) + qJ(5);
	t13 = cos(t15);
	t12 = sin(t15);
	t29 = t12 * t19;
	t5 = t10 * t29 + t11 * t13;
	t28 = t13 * t19;
	t6 = -t10 * t28 + t11 * t12;
	t33 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t7 = t10 * t13 - t11 * t29;
	t8 = t10 * t12 + t11 * t28;
	t32 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t16 = sin(qJ(4));
	t31 = pkin(4) * t16;
	t27 = t16 * t19;
	t25 = pkin(7) + t31;
	t23 = -r_i_i_C(1) * t12 - r_i_i_C(2) * t13;
	t22 = r_i_i_C(1) * t13 - r_i_i_C(2) * t12 + t9;
	t21 = -t22 * t17 + t30 * t19;
	t1 = [-sin(qJ(1)) * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t25 * t11 - t34 * t10, 0, t21 * t11, (t10 * t18 - t11 * t27) * pkin(4) + t32, t32, 0; cos(qJ(1)) * pkin(1) + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t25 * t10 + t34 * t11, 0, t21 * t10, (-t10 * t27 - t11 * t18) * pkin(4) + t33, t33, 0; 0, 1, t22 * t19 + t24, (t23 - t31) * t17, t23 * t17, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (252->34), mult. (205->47), div. (0->0), fcn. (234->10), ass. (0->32)
	t42 = pkin(5) + r_i_i_C(1);
	t35 = r_i_i_C(3) + qJ(6);
	t24 = cos(qJ(4));
	t15 = t24 * pkin(4) + pkin(3);
	t25 = cos(qJ(3));
	t23 = sin(qJ(3));
	t40 = r_i_i_C(2) + pkin(9) + pkin(8);
	t31 = t40 * t23;
	t44 = t25 * t15 + pkin(2) + t31;
	t21 = qJ(4) + qJ(5);
	t18 = sin(t21);
	t19 = cos(t21);
	t43 = t35 * t18 + t42 * t19 + t15;
	t22 = sin(qJ(4));
	t41 = pkin(4) * t22;
	t39 = t18 * t25;
	t37 = t19 * t25;
	t36 = t22 * t25;
	t34 = t35 * t19 * t23;
	t33 = t42 * t18;
	t32 = pkin(7) + t41;
	t20 = qJ(1) + pkin(10);
	t16 = sin(t20);
	t17 = cos(t20);
	t7 = t16 * t39 + t17 * t19;
	t8 = t16 * t37 - t17 * t18;
	t29 = t35 * t8 - t42 * t7;
	t10 = t16 * t18 + t17 * t37;
	t9 = -t16 * t19 + t17 * t39;
	t28 = t35 * t10 - t42 * t9;
	t27 = -t43 * t23 + t40 * t25;
	t1 = [-sin(qJ(1)) * pkin(1) - t42 * t8 - t35 * t7 + t32 * t17 - t44 * t16, 0, t27 * t17, (t16 * t24 - t17 * t36) * pkin(4) + t28, t28, t9; cos(qJ(1)) * pkin(1) + t35 * t9 + t32 * t16 + t42 * t10 + t44 * t17, 0, t27 * t16, (-t16 * t36 - t17 * t24) * pkin(4) + t29, t29, t7; 0, 1, t43 * t25 + t31, (-t33 - t41) * t23 + t34, -t23 * t33 + t34, t23 * t18;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end