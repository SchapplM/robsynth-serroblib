% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRP7
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
% Datum: 2019-10-10 10:39
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:39:16
	% EndTime: 2019-10-10 10:39:16
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:39:16
	% EndTime: 2019-10-10 10:39:16
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
	% StartTime: 2019-10-10 10:39:16
	% EndTime: 2019-10-10 10:39:16
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
	% StartTime: 2019-10-10 10:39:16
	% EndTime: 2019-10-10 10:39:16
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (22->9), mult. (44->12), div. (0->0), fcn. (47->4), ass. (0->11)
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(3) + qJ(3);
	t9 = pkin(2) + r_i_i_C(1);
	t6 = t7 * t1 + t9 * t3;
	t10 = pkin(1) + t6;
	t8 = pkin(7) + r_i_i_C(2);
	t5 = -t9 * t1 + t7 * t3;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t11 = [-t10 * t2 + t8 * t4, t5 * t4, t4 * t1, 0, 0, 0; t10 * t4 + t8 * t2, t5 * t2, t2 * t1, 0, 0, 0; 0, t6, -t3, 0, 0, 0;];
	Ja_transl = t11;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:39:16
	% EndTime: 2019-10-10 10:39:16
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (48->21), mult. (109->31), div. (0->0), fcn. (128->6), ass. (0->23)
	t12 = cos(qJ(2));
	t24 = pkin(2) + pkin(3);
	t9 = sin(qJ(2));
	t14 = t9 * qJ(3) + t24 * t12;
	t25 = pkin(1) + t14;
	t8 = sin(qJ(4));
	t23 = t12 * t8;
	t13 = cos(qJ(1));
	t22 = t13 * t9;
	t20 = pkin(7) - pkin(8) - r_i_i_C(3);
	t10 = sin(qJ(1));
	t11 = cos(qJ(4));
	t6 = t11 * t9 - t23;
	t1 = t6 * t10;
	t16 = t11 * t12 + t8 * t9;
	t2 = t16 * t10;
	t19 = t1 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t3 = -t11 * t22 + t13 * t23;
	t4 = t16 * t13;
	t18 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t17 = -r_i_i_C(1) * t16 - t6 * r_i_i_C(2);
	t15 = qJ(3) * t12 - t24 * t9;
	t5 = [-t2 * r_i_i_C(1) - t1 * r_i_i_C(2) - t25 * t10 + t20 * t13, t15 * t13 - t18, t22, t18, 0, 0; t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t20 * t10 + t25 * t13, t15 * t10 - t19, t10 * t9, t19, 0, 0; 0, t14 - t17, -t12, t17, 0, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:39:16
	% EndTime: 2019-10-10 10:39:16
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (111->35), mult. (271->53), div. (0->0), fcn. (328->8), ass. (0->30)
	t15 = sin(qJ(2));
	t18 = cos(qJ(2));
	t33 = pkin(2) + pkin(3);
	t20 = t15 * qJ(3) + t33 * t18;
	t37 = pkin(1) + t20;
	t29 = sin(qJ(4));
	t30 = cos(qJ(4));
	t8 = t15 * t30 - t18 * t29;
	t14 = sin(qJ(5));
	t17 = cos(qJ(5));
	t22 = t17 * r_i_i_C(1) - t14 * r_i_i_C(2) + pkin(4);
	t16 = sin(qJ(1));
	t3 = t8 * t16;
	t31 = -r_i_i_C(3) - pkin(9);
	t7 = t15 * t29 + t18 * t30;
	t4 = t7 * t16;
	t36 = t22 * t3 - t31 * t4;
	t19 = cos(qJ(1));
	t24 = t19 * t29;
	t25 = t19 * t30;
	t5 = -t15 * t24 - t18 * t25;
	t6 = -t15 * t25 + t18 * t24;
	t35 = -t22 * t6 + t31 * t5;
	t34 = -t22 * t7 - t31 * t8;
	t32 = pkin(7) - pkin(8);
	t23 = -t14 * r_i_i_C(1) - t17 * r_i_i_C(2);
	t21 = qJ(3) * t18 - t33 * t15;
	t2 = -t16 * t14 - t5 * t17;
	t1 = t5 * t14 - t16 * t17;
	t9 = [-t31 * t3 - t22 * t4 + (t23 + t32) * t19 - t37 * t16, t21 * t19 - t35, t19 * t15, t35, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; -t5 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t32 * t16 + t37 * t19 - t31 * t6, t21 * t16 - t36, t16 * t15, t36, (-t4 * t14 + t19 * t17) * r_i_i_C(1) + (-t19 * t14 - t4 * t17) * r_i_i_C(2), 0; 0, t20 - t34, -t18, t34, t23 * t8, 0;];
	Ja_transl = t9;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:39:16
	% EndTime: 2019-10-10 10:39:16
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (168->36), mult. (414->52), div. (0->0), fcn. (509->8), ass. (0->31)
	t19 = sin(qJ(5));
	t21 = cos(qJ(5));
	t34 = r_i_i_C(3) + qJ(6);
	t40 = pkin(5) + r_i_i_C(1);
	t25 = t34 * t19 + t40 * t21 + pkin(4);
	t38 = pkin(9) + r_i_i_C(2);
	t22 = cos(qJ(2));
	t35 = sin(qJ(4));
	t36 = sin(qJ(2));
	t37 = cos(qJ(4));
	t12 = -t22 * t35 + t36 * t37;
	t20 = sin(qJ(1));
	t7 = t12 * t20;
	t11 = t22 * t37 + t36 * t35;
	t8 = t11 * t20;
	t48 = t25 * t7 + t38 * t8;
	t23 = cos(qJ(1));
	t10 = t12 * t23;
	t9 = t11 * t23;
	t47 = t25 * t10 + t38 * t9;
	t46 = -t25 * t11 + t38 * t12;
	t41 = pkin(2) + pkin(3);
	t27 = t36 * qJ(3) + t41 * t22;
	t43 = pkin(1) + t27;
	t39 = pkin(7) - pkin(8);
	t1 = t8 * t19 - t23 * t21;
	t28 = t23 * t19 + t8 * t21;
	t24 = qJ(3) * t22 - t36 * t41;
	t6 = -t20 * t19 + t9 * t21;
	t5 = t9 * t19 + t20 * t21;
	t2 = [-t8 * pkin(4) - t34 * t1 - t20 * t43 + t39 * t23 - t40 * t28 + t38 * t7, t24 * t23 - t47, t23 * t36, t47, t34 * t6 - t40 * t5, t5; t9 * pkin(4) - t38 * t10 + t39 * t20 + t23 * t43 + t34 * t5 + t40 * t6, t24 * t20 - t48, t20 * t36, t48, -t40 * t1 + t34 * t28, t1; 0, t27 - t46, -t22, t46, (-t40 * t19 + t34 * t21) * t12, t12 * t19;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end