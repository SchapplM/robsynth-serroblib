% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
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
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
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
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
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
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (104->31), mult. (255->53), div. (0->0), fcn. (319->8), ass. (0->27)
	t18 = sin(qJ(3));
	t21 = cos(qJ(3));
	t29 = r_i_i_C(3) + qJ(4);
	t34 = pkin(3) + r_i_i_C(1);
	t35 = t29 * t18 + t34 * t21 + pkin(2);
	t33 = pkin(9) + r_i_i_C(2);
	t17 = sin(pkin(6));
	t20 = sin(qJ(1));
	t32 = t17 * t20;
	t31 = t17 * t21;
	t23 = cos(qJ(1));
	t30 = t17 * t23;
	t28 = cos(pkin(6));
	t19 = sin(qJ(2));
	t22 = cos(qJ(2));
	t25 = t23 * t28;
	t10 = t19 * t25 + t20 * t22;
	t27 = t10 * t21 - t18 * t30;
	t26 = t20 * t28;
	t1 = t10 * t18 + t21 * t30;
	t12 = -t19 * t26 + t23 * t22;
	t11 = t23 * t19 + t22 * t26;
	t9 = t20 * t19 - t22 * t25;
	t7 = t17 * t19 * t18 - t28 * t21;
	t6 = t12 * t21 + t18 * t32;
	t5 = t12 * t18 - t20 * t31;
	t2 = [-t20 * pkin(1) - t10 * pkin(2) + pkin(8) * t30 - t29 * t1 - t34 * t27 - t33 * t9, -t11 * t35 + t33 * t12, t29 * t6 - t34 * t5, t5, 0, 0; t23 * pkin(1) + t12 * pkin(2) + pkin(8) * t32 + t33 * t11 + t29 * t5 + t34 * t6, t33 * t10 - t35 * t9, -t34 * t1 + t29 * t27, t1, 0, 0; 0, (t33 * t19 + t35 * t22) * t17, t29 * (t28 * t18 + t19 * t31) - t34 * t7, t7, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (134->33), mult. (325->54), div. (0->0), fcn. (409->8), ass. (0->27)
	t18 = sin(qJ(3));
	t21 = cos(qJ(3));
	t29 = pkin(3) + pkin(4) - r_i_i_C(2);
	t31 = r_i_i_C(1) + qJ(4);
	t35 = t31 * t18 + t29 * t21 + pkin(2);
	t17 = sin(pkin(6));
	t20 = sin(qJ(1));
	t34 = t17 * t20;
	t33 = t17 * t21;
	t23 = cos(qJ(1));
	t32 = t17 * t23;
	t30 = cos(pkin(6));
	t28 = pkin(9) - r_i_i_C(3) - qJ(5);
	t27 = t20 * t30;
	t26 = t23 * t30;
	t19 = sin(qJ(2));
	t22 = cos(qJ(2));
	t10 = t19 * t26 + t20 * t22;
	t1 = t10 * t18 + t21 * t32;
	t25 = -t10 * t21 + t18 * t32;
	t12 = -t19 * t27 + t23 * t22;
	t11 = t23 * t19 + t22 * t27;
	t9 = t20 * t19 - t22 * t26;
	t7 = t17 * t19 * t18 - t30 * t21;
	t6 = t12 * t21 + t18 * t34;
	t5 = t12 * t18 - t20 * t33;
	t2 = [-t20 * pkin(1) - t10 * pkin(2) + pkin(8) * t32 - t31 * t1 + t29 * t25 - t28 * t9, -t11 * t35 + t28 * t12, -t29 * t5 + t31 * t6, t5, -t11, 0; t23 * pkin(1) + t12 * pkin(2) + pkin(8) * t34 + t28 * t11 + t29 * t6 + t31 * t5, t28 * t10 - t35 * t9, -t29 * t1 - t31 * t25, t1, -t9, 0; 0, (t28 * t19 + t35 * t22) * t17, t31 * (t30 * t18 + t19 * t33) - t29 * t7, t7, t17 * t22, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (223->48), mult. (556->78), div. (0->0), fcn. (707->10), ass. (0->34)
	t21 = sin(qJ(3));
	t25 = cos(qJ(3));
	t20 = sin(qJ(6));
	t24 = cos(qJ(6));
	t37 = pkin(5) + qJ(4);
	t30 = t24 * r_i_i_C(1) - t20 * r_i_i_C(2) + t37;
	t34 = pkin(3) + pkin(4) + pkin(10) + r_i_i_C(3);
	t42 = t30 * t21 + t34 * t25 + pkin(2);
	t41 = cos(qJ(1));
	t19 = sin(pkin(6));
	t23 = sin(qJ(1));
	t40 = t19 * t23;
	t39 = t19 * t25;
	t26 = cos(qJ(2));
	t38 = t19 * t26;
	t36 = pkin(9) - qJ(5);
	t35 = cos(pkin(6));
	t33 = t19 * t41;
	t32 = t23 * t35;
	t31 = t35 * t41;
	t29 = t20 * r_i_i_C(1) + t24 * r_i_i_C(2) - t36;
	t22 = sin(qJ(2));
	t12 = t22 * t31 + t23 * t26;
	t3 = t12 * t21 + t25 * t33;
	t28 = t12 * t25 - t21 * t33;
	t14 = -t22 * t32 + t41 * t26;
	t13 = t41 * t22 + t26 * t32;
	t11 = t23 * t22 - t26 * t31;
	t9 = t19 * t22 * t21 - t35 * t25;
	t8 = t14 * t25 + t21 * t40;
	t7 = t14 * t21 - t23 * t39;
	t2 = -t13 * t20 + t7 * t24;
	t1 = -t13 * t24 - t7 * t20;
	t4 = [-t23 * pkin(1) - t12 * pkin(2) + pkin(8) * t33 + t29 * t11 - t34 * t28 - t30 * t3, -t13 * t42 - t29 * t14, t30 * t8 - t34 * t7, t7, -t13, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t41 * pkin(1) + t14 * pkin(2) + pkin(8) * t40 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t36 * t13 + t34 * t8 + t37 * t7, -t11 * t42 - t29 * t12, t30 * t28 - t34 * t3, t3, -t11, (-t11 * t24 - t3 * t20) * r_i_i_C(1) + (t11 * t20 - t3 * t24) * r_i_i_C(2); 0, (-t29 * t22 + t42 * t26) * t19, -t34 * t9 + t30 * (t35 * t21 + t22 * t39), t9, t38, (-t9 * t20 + t24 * t38) * r_i_i_C(1) + (-t20 * t38 - t9 * t24) * r_i_i_C(2);];
	Ja_transl = t4;
else
	Ja_transl=NaN(3,6);
end