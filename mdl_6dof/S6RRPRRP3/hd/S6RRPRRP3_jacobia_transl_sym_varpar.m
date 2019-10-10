% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
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
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
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
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (74->25), mult. (90->34), div. (0->0), fcn. (100->8), ass. (0->22)
	t24 = pkin(8) + r_i_i_C(3);
	t9 = qJ(2) + pkin(10);
	t6 = sin(t9);
	t26 = t24 * t6 + cos(qJ(2)) * pkin(2);
	t7 = cos(t9);
	t25 = t7 * pkin(3) + pkin(1) + t26;
	t11 = sin(qJ(4));
	t15 = cos(qJ(1));
	t23 = t11 * t15;
	t13 = sin(qJ(1));
	t22 = t13 * t11;
	t14 = cos(qJ(4));
	t21 = t13 * t14;
	t20 = t14 * t15;
	t17 = r_i_i_C(1) * t14 - r_i_i_C(2) * t11 + pkin(3);
	t16 = -sin(qJ(2)) * pkin(2) - t17 * t6 + t24 * t7;
	t10 = -qJ(3) - pkin(7);
	t4 = t7 * t20 + t22;
	t3 = -t7 * t23 + t21;
	t2 = -t7 * t21 + t23;
	t1 = t7 * t22 + t20;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t15 * t10 - t25 * t13, t16 * t15, t13, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t13 * t10 + t25 * t15, t16 * t13, -t15, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, t17 * t7 + t26, 0, (-r_i_i_C(1) * t11 - r_i_i_C(2) * t14) * t6, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (141->32), mult. (131->44), div. (0->0), fcn. (145->10), ass. (0->31)
	t16 = qJ(2) + pkin(10);
	t11 = sin(t16);
	t36 = r_i_i_C(3) + pkin(9) + pkin(8);
	t41 = cos(qJ(2)) * pkin(2) + t36 * t11;
	t12 = cos(t16);
	t22 = cos(qJ(4));
	t9 = t22 * pkin(4) + pkin(3);
	t40 = t12 * t9 + pkin(1) + t41;
	t17 = qJ(4) + qJ(5);
	t14 = cos(t17);
	t23 = cos(qJ(1));
	t31 = t23 * t14;
	t13 = sin(t17);
	t21 = sin(qJ(1));
	t34 = t21 * t13;
	t5 = t12 * t34 + t31;
	t32 = t23 * t13;
	t33 = t21 * t14;
	t6 = -t12 * t33 + t32;
	t39 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t7 = -t12 * t32 + t33;
	t8 = t12 * t31 + t34;
	t38 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t19 = sin(qJ(4));
	t37 = pkin(4) * t19;
	t35 = t12 * t19;
	t29 = qJ(3) + pkin(7) + t37;
	t27 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
	t26 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
	t25 = -sin(qJ(2)) * pkin(2) - t26 * t11 + t36 * t12;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t40 * t21 + t29 * t23, t25 * t23, t21, (t21 * t22 - t23 * t35) * pkin(4) + t38, t38, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t29 * t21 + t40 * t23, t25 * t21, -t23, (-t21 * t35 - t22 * t23) * pkin(4) + t39, t39, 0; 0, t12 * t26 + t41, 0, (t27 - t37) * t11, t27 * t11, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:33:36
	% EndTime: 2019-10-10 10:33:36
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (187->39), mult. (157->49), div. (0->0), fcn. (174->10), ass. (0->30)
	t21 = qJ(2) + pkin(10);
	t14 = sin(t21);
	t37 = r_i_i_C(3) + qJ(6) + pkin(9) + pkin(8);
	t42 = cos(qJ(2)) * pkin(2) + t37 * t14;
	t15 = cos(t21);
	t22 = qJ(4) + qJ(5);
	t17 = cos(t22);
	t11 = pkin(5) * t17 + cos(qJ(4)) * pkin(4);
	t9 = pkin(3) + t11;
	t41 = t15 * t9 + pkin(1) + t42;
	t16 = sin(t22);
	t25 = sin(qJ(1));
	t33 = t25 * t16;
	t26 = cos(qJ(1));
	t34 = t17 * t26;
	t5 = t15 * t33 + t34;
	t32 = t25 * t17;
	t35 = t16 * t26;
	t6 = -t15 * t32 + t35;
	t40 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t7 = -t15 * t35 + t32;
	t8 = t15 * t34 + t33;
	t39 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t38 = r_i_i_C(2) * t17;
	t10 = pkin(5) * t16 + sin(qJ(4)) * pkin(4);
	t36 = t10 * t15;
	t31 = t10 + qJ(3) + pkin(7);
	t28 = r_i_i_C(1) * t17 - r_i_i_C(2) * t16 + t9;
	t27 = -sin(qJ(2)) * pkin(2) - t28 * t14 + t37 * t15;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t41 * t25 + t31 * t26, t27 * t26, t25, t25 * t11 - t26 * t36 + t39, t7 * pkin(5) + t39, t26 * t14; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t31 * t25 + t41 * t26, t27 * t25, -t26, -t11 * t26 - t25 * t36 + t40, -t5 * pkin(5) + t40, t25 * t14; 0, t28 * t15 + t42, 0, (-r_i_i_C(1) * t16 - t10 - t38) * t14, (-t38 + (-pkin(5) - r_i_i_C(1)) * t16) * t14, -t15;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end