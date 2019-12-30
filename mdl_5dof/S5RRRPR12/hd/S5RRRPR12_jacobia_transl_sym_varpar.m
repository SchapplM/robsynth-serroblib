% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR12
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRPR12_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR12_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR12_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:19:07
	% EndTime: 2019-12-29 20:19:07
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:19:13
	% EndTime: 2019-12-29 20:19:13
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:19:07
	% EndTime: 2019-12-29 20:19:08
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (19->14), mult. (46->24), div. (0->0), fcn. (54->6), ass. (0->16)
	t7 = sin(qJ(2));
	t8 = sin(qJ(1));
	t15 = t8 * t7;
	t9 = cos(qJ(2));
	t14 = t8 * t9;
	t10 = cos(qJ(1));
	t13 = t10 * t7;
	t12 = t10 * t9;
	t5 = sin(pkin(5));
	t11 = (pkin(7) + r_i_i_C(3)) * t5;
	t6 = cos(pkin(5));
	t4 = -t15 * t6 + t12;
	t3 = -t14 * t6 - t13;
	t2 = -t13 * t6 - t14;
	t1 = -t12 * t6 + t15;
	t16 = [-t8 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t10 * t11, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0; t10 * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t11 * t8, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0; 0, (r_i_i_C(1) * t9 - r_i_i_C(2) * t7) * t5, 0, 0, 0;];
	Ja_transl = t16;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:19:07
	% EndTime: 2019-12-29 20:19:08
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (61->29), mult. (150->51), div. (0->0), fcn. (184->8), ass. (0->26)
	t27 = -r_i_i_C(3) - pkin(8);
	t13 = sin(qJ(1));
	t9 = sin(pkin(5));
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
	t10 = cos(pkin(5));
	t4 = t10 * t22 + t21;
	t19 = t11 * t24 - t14 * t4;
	t18 = t14 * r_i_i_C(1) - t11 * r_i_i_C(2) + pkin(2);
	t17 = t4 * t11 + t14 * t24;
	t6 = -t10 * t23 + t20;
	t5 = t10 * t21 + t22;
	t3 = -t10 * t20 + t23;
	t2 = t11 * t26 + t14 * t6;
	t1 = -t11 * t6 + t13 * t25;
	t7 = [-t13 * pkin(1) - t4 * pkin(2) + pkin(7) * t24 + t19 * r_i_i_C(1) + t17 * r_i_i_C(2) + t27 * t3, -t18 * t5 - t27 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0; pkin(1) * t16 + t6 * pkin(2) + pkin(7) * t26 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t27 * t5, -t18 * t3 - t27 * t4, -t17 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0; 0, (-t27 * t12 + t18 * t15) * t9, (-t11 * t12 * t9 + t10 * t14) * r_i_i_C(1) + (-t10 * t11 - t12 * t25) * r_i_i_C(2), 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:19:12
	% EndTime: 2019-12-29 20:19:13
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (132->41), mult. (334->69), div. (0->0), fcn. (423->10), ass. (0->29)
	t20 = sin(qJ(3));
	t23 = cos(qJ(3));
	t17 = sin(pkin(10));
	t19 = cos(pkin(10));
	t28 = r_i_i_C(1) * t19 - r_i_i_C(2) * t17 + pkin(3);
	t33 = r_i_i_C(3) + qJ(4);
	t37 = t33 * t20 + t28 * t23 + pkin(2);
	t18 = sin(pkin(5));
	t22 = sin(qJ(1));
	t36 = t18 * t22;
	t35 = t18 * t23;
	t25 = cos(qJ(1));
	t34 = t18 * t25;
	t32 = cos(pkin(5));
	t21 = sin(qJ(2));
	t24 = cos(qJ(2));
	t29 = t25 * t32;
	t10 = t21 * t29 + t22 * t24;
	t31 = t10 * t23 - t20 * t34;
	t30 = t22 * t32;
	t27 = t17 * r_i_i_C(1) + t19 * r_i_i_C(2) + pkin(8);
	t1 = t10 * t20 + t23 * t34;
	t12 = -t21 * t30 + t25 * t24;
	t11 = t25 * t21 + t24 * t30;
	t9 = t22 * t21 - t24 * t29;
	t7 = t18 * t21 * t20 - t32 * t23;
	t6 = t12 * t23 + t20 * t36;
	t5 = t12 * t20 - t22 * t35;
	t2 = [(-t9 * t17 - t19 * t31) * r_i_i_C(1) + (t17 * t31 - t9 * t19) * r_i_i_C(2) - t31 * pkin(3) - t10 * pkin(2) - t9 * pkin(8) - t22 * pkin(1) + pkin(7) * t34 - t33 * t1, -t11 * t37 + t27 * t12, -t28 * t5 + t33 * t6, t5, 0; (t11 * t17 + t6 * t19) * r_i_i_C(1) + (t11 * t19 - t6 * t17) * r_i_i_C(2) + t6 * pkin(3) + t12 * pkin(2) + t11 * pkin(8) + t25 * pkin(1) + pkin(7) * t36 + t33 * t5, t27 * t10 - t37 * t9, -t28 * t1 + t33 * t31, t1, 0; 0, (t27 * t21 + t37 * t24) * t18, t33 * (t32 * t20 + t21 * t35) - t28 * t7, t7, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:19:07
	% EndTime: 2019-12-29 20:19:08
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (218->47), mult. (427->80), div. (0->0), fcn. (540->12), ass. (0->36)
	t26 = sin(qJ(3));
	t29 = cos(qJ(3));
	t19 = cos(pkin(10)) * pkin(4) + pkin(3);
	t22 = pkin(10) + qJ(5);
	t20 = sin(t22);
	t21 = cos(t22);
	t33 = t21 * r_i_i_C(1) - t20 * r_i_i_C(2) + t19;
	t43 = r_i_i_C(3) + pkin(9) + qJ(4);
	t44 = t43 * t26 + t33 * t29 + pkin(2);
	t42 = cos(qJ(1));
	t24 = sin(pkin(5));
	t28 = sin(qJ(1));
	t41 = t24 * t28;
	t40 = t24 * t29;
	t30 = cos(qJ(2));
	t39 = t24 * t30;
	t38 = cos(pkin(5));
	t37 = sin(pkin(10)) * pkin(4) + pkin(8);
	t36 = t24 * t42;
	t27 = sin(qJ(2));
	t34 = t38 * t42;
	t12 = t27 * t34 + t28 * t30;
	t4 = t12 * t29 - t26 * t36;
	t35 = t28 * t38;
	t3 = t12 * t26 + t29 * t36;
	t32 = t20 * r_i_i_C(1) + t21 * r_i_i_C(2) + t37;
	t14 = -t27 * t35 + t42 * t30;
	t13 = t42 * t27 + t30 * t35;
	t11 = t28 * t27 - t30 * t34;
	t10 = t38 * t26 + t27 * t40;
	t9 = t24 * t27 * t26 - t38 * t29;
	t8 = t14 * t29 + t26 * t41;
	t7 = t14 * t26 - t28 * t40;
	t2 = t13 * t20 + t8 * t21;
	t1 = t13 * t21 - t8 * t20;
	t5 = [-t28 * pkin(1) - t12 * pkin(2) + pkin(7) * t36 - t32 * t11 - t43 * t3 - t33 * t4, -t13 * t44 + t32 * t14, -t33 * t7 + t43 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t42 * pkin(1) + t14 * pkin(2) + pkin(7) * t41 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t37 * t13 + t8 * t19 + t43 * t7, -t11 * t44 + t32 * t12, -t33 * t3 + t43 * t4, t3, (t11 * t21 - t4 * t20) * r_i_i_C(1) + (-t11 * t20 - t4 * t21) * r_i_i_C(2); 0, (t32 * t27 + t44 * t30) * t24, t43 * t10 - t33 * t9, t9, (-t10 * t20 - t21 * t39) * r_i_i_C(1) + (-t10 * t21 + t20 * t39) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end