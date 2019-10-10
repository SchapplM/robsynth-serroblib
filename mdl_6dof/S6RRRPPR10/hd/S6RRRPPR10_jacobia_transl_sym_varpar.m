% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR10_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR10_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR10_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (19->14), mult. (46->24), div. (0->0), fcn. (54->6), ass. (0->16)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t15 = t8 * t7;
	t9 = cos(qJ(2));
	t14 = t8 * t9;
	t10 = cos(qJ(1));
	t13 = t10 * t7;
	t12 = t10 * t9;
	t5 = sin(pkin(6));
	t11 = (pkin(8) + r_i_i_C(3)) * t5;
	t6 = cos(pkin(6));
	t4 = -t6 * t15 + t12;
	t3 = -t6 * t14 - t13;
	t2 = -t6 * t13 - t14;
	t1 = -t6 * t12 + t15;
	t16 = [-t8 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t10 * t11, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0, 0; t10 * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t11 * t8, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0; 0, (r_i_i_C(1) * t9 - r_i_i_C(2) * t7) * t5, 0, 0, 0, 0;];
	Ja_transl = t16;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (61->29), mult. (150->51), div. (0->0), fcn. (184->8), ass. (0->26)
	t27 = -r_i_i_C(3) - pkin(9);
	t13 = sin(qJ(1));
	t9 = sin(pkin(6));
	t26 = t13 * t9;
	t14 = cos(qJ(3));
	t25 = t14 * t9;
	t16 = cos(qJ(1));
	t24 = t16 * t9;
	t12 = sin(qJ(2));
	t23 = t12 * t13;
	t22 = t12 * t16;
	t15 = cos(qJ(2));
	t21 = t13 * t15;
	t20 = t15 * t16;
	t11 = sin(qJ(3));
	t10 = cos(pkin(6));
	t4 = t10 * t22 + t21;
	t19 = t11 * t24 - t14 * t4;
	t18 = t14 * r_i_i_C(1) - t11 * r_i_i_C(2) + pkin(2);
	t17 = t4 * t11 + t14 * t24;
	t6 = -t10 * t23 + t20;
	t5 = t10 * t21 + t22;
	t3 = -t10 * t20 + t23;
	t2 = t11 * t26 + t14 * t6;
	t1 = -t11 * t6 + t13 * t25;
	t7 = [-t13 * pkin(1) - t4 * pkin(2) + pkin(8) * t24 + t19 * r_i_i_C(1) + t17 * r_i_i_C(2) + t27 * t3, -t18 * t5 - t27 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; pkin(1) * t16 + t6 * pkin(2) + pkin(8) * t26 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t27 * t5, -t18 * t3 - t27 * t4, -t17 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, (-t27 * t12 + t18 * t15) * t9, (-t11 * t12 * t9 + t10 * t14) * r_i_i_C(1) + (-t10 * t11 - t12 * t25) * r_i_i_C(2), 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (104->31), mult. (255->53), div. (0->0), fcn. (319->8), ass. (0->27)
	t16 = sin(qJ(3));
	t19 = cos(qJ(3));
	t27 = r_i_i_C(3) + qJ(4);
	t32 = pkin(3) - r_i_i_C(2);
	t33 = t27 * t16 + t32 * t19 + pkin(2);
	t31 = pkin(9) + r_i_i_C(1);
	t15 = sin(pkin(6));
	t18 = sin(qJ(1));
	t30 = t15 * t18;
	t29 = t15 * t19;
	t21 = cos(qJ(1));
	t28 = t15 * t21;
	t26 = cos(pkin(6));
	t25 = t18 * t26;
	t24 = t21 * t26;
	t17 = sin(qJ(2));
	t20 = cos(qJ(2));
	t10 = t17 * t24 + t18 * t20;
	t1 = t10 * t16 + t19 * t28;
	t23 = -t10 * t19 + t16 * t28;
	t12 = -t17 * t25 + t21 * t20;
	t11 = t21 * t17 + t20 * t25;
	t9 = t18 * t17 - t20 * t24;
	t7 = t15 * t17 * t16 - t26 * t19;
	t6 = t12 * t19 + t16 * t30;
	t5 = t12 * t16 - t18 * t29;
	t2 = [-t18 * pkin(1) - t10 * pkin(2) + pkin(8) * t28 - t27 * t1 + t32 * t23 - t31 * t9, -t11 * t33 + t31 * t12, t27 * t6 - t32 * t5, t5, 0, 0; t21 * pkin(1) + t12 * pkin(2) + pkin(8) * t30 + t31 * t11 + t27 * t5 + t32 * t6, t31 * t10 - t33 * t9, -t32 * t1 - t27 * t23, t1, 0, 0; 0, (t31 * t17 + t33 * t20) * t15, t27 * (t26 * t16 + t17 * t29) - t32 * t7, t7, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:19
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (165->34), mult. (412->57), div. (0->0), fcn. (524->10), ass. (0->30)
	t19 = sin(qJ(3));
	t22 = cos(qJ(3));
	t16 = sin(pkin(11));
	t18 = cos(pkin(11));
	t26 = r_i_i_C(1) * t16 + r_i_i_C(2) * t18 + qJ(4);
	t30 = pkin(3) + r_i_i_C(3) + qJ(5);
	t35 = t26 * t19 + t30 * t22 + pkin(2);
	t34 = cos(qJ(1));
	t17 = sin(pkin(6));
	t21 = sin(qJ(1));
	t33 = t17 * t21;
	t32 = t17 * t22;
	t31 = cos(pkin(6));
	t29 = t17 * t34;
	t28 = t21 * t31;
	t27 = t31 * t34;
	t25 = t18 * r_i_i_C(1) - t16 * r_i_i_C(2) + pkin(4) + pkin(9);
	t20 = sin(qJ(2));
	t23 = cos(qJ(2));
	t10 = t20 * t27 + t21 * t23;
	t1 = t10 * t19 + t22 * t29;
	t2 = t10 * t22 - t19 * t29;
	t12 = -t20 * t28 + t34 * t23;
	t11 = t34 * t20 + t23 * t28;
	t9 = t21 * t20 - t23 * t27;
	t8 = t31 * t19 + t20 * t32;
	t7 = t17 * t20 * t19 - t31 * t22;
	t6 = t12 * t22 + t19 * t33;
	t5 = t12 * t19 - t21 * t32;
	t3 = [-t21 * pkin(1) - t10 * pkin(2) + pkin(8) * t29 - t26 * t1 - t30 * t2 - t25 * t9, -t11 * t35 + t25 * t12, t26 * t6 - t30 * t5, t5, t6, 0; t34 * pkin(1) + t12 * pkin(2) + pkin(8) * t33 + t25 * t11 + t26 * t5 + t30 * t6, t25 * t10 - t35 * t9, -t30 * t1 + t26 * t2, t1, t2, 0; 0, (t25 * t20 + t35 * t23) * t17, t26 * t8 - t30 * t7, t7, t8, 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:33:19
	% EndTime: 2019-10-10 11:33:20
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (258->48), mult. (532->80), div. (0->0), fcn. (675->12), ass. (0->36)
	t25 = sin(qJ(3));
	t28 = cos(qJ(3));
	t21 = pkin(11) + qJ(6);
	t19 = sin(t21);
	t20 = cos(t21);
	t35 = pkin(5) * sin(pkin(11)) + qJ(4);
	t31 = t19 * r_i_i_C(1) + t20 * r_i_i_C(2) + t35;
	t37 = pkin(3) + r_i_i_C(3) + pkin(10) + qJ(5);
	t44 = t25 * t31 + t28 * t37 + pkin(2);
	t43 = pkin(9) + cos(pkin(11)) * pkin(5) + pkin(4);
	t42 = cos(qJ(1));
	t23 = sin(pkin(6));
	t27 = sin(qJ(1));
	t41 = t23 * t27;
	t40 = t23 * t28;
	t29 = cos(qJ(2));
	t39 = t23 * t29;
	t38 = cos(pkin(6));
	t36 = t23 * t42;
	t34 = t27 * t38;
	t33 = t38 * t42;
	t32 = t20 * r_i_i_C(1) - t19 * r_i_i_C(2) + t43;
	t26 = sin(qJ(2));
	t12 = t26 * t33 + t27 * t29;
	t3 = t12 * t25 + t28 * t36;
	t4 = t12 * t28 - t25 * t36;
	t14 = -t26 * t34 + t29 * t42;
	t13 = t26 * t42 + t29 * t34;
	t11 = t27 * t26 - t29 * t33;
	t10 = t25 * t38 + t26 * t40;
	t9 = t23 * t26 * t25 - t28 * t38;
	t8 = t14 * t28 + t25 * t41;
	t7 = t14 * t25 - t27 * t40;
	t2 = t13 * t20 + t7 * t19;
	t1 = -t13 * t19 + t7 * t20;
	t5 = [-t27 * pkin(1) - t12 * pkin(2) + pkin(8) * t36 - t11 * t32 - t3 * t31 - t37 * t4, -t13 * t44 + t14 * t32, t31 * t8 - t37 * t7, t7, t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t42 * pkin(1) + t14 * pkin(2) + pkin(8) * t41 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t43 * t13 + t35 * t7 + t37 * t8, -t11 * t44 + t12 * t32, -t3 * t37 + t31 * t4, t3, t4, (-t11 * t19 + t3 * t20) * r_i_i_C(1) + (-t11 * t20 - t3 * t19) * r_i_i_C(2); 0, (t32 * t26 + t44 * t29) * t23, t10 * t31 - t37 * t9, t9, t10, (t19 * t39 + t9 * t20) * r_i_i_C(1) + (-t9 * t19 + t20 * t39) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end