% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR3
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:39
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
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
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (49->24), mult. (114->41), div. (0->0), fcn. (143->8), ass. (0->19)
	t16 = cos(qJ(2));
	t10 = sin(pkin(11));
	t12 = cos(pkin(11));
	t14 = sin(qJ(2));
	t6 = t10 * t14 - t16 * t12;
	t24 = -t16 * pkin(2) + t6 * r_i_i_C(1);
	t13 = cos(pkin(6));
	t21 = t13 * t16;
	t20 = -pkin(1) + t24;
	t19 = t10 * t16 + t12 * t14;
	t11 = sin(pkin(6));
	t4 = t19 * t13;
	t18 = -pkin(2) * t13 * t14 - t4 * r_i_i_C(1) + (r_i_i_C(3) + pkin(8) + qJ(3)) * t11;
	t17 = cos(qJ(1));
	t15 = sin(qJ(1));
	t3 = t6 * t13;
	t2 = t15 * t3 - t17 * t19;
	t1 = -t15 * t19 - t17 * t3;
	t5 = [-t1 * r_i_i_C(2) + t20 * t15 + t18 * t17, t2 * r_i_i_C(1) + (t15 * t4 + t17 * t6) * r_i_i_C(2) + (-t14 * t17 - t15 * t21) * pkin(2), t15 * t11, 0, 0, 0; t2 * r_i_i_C(2) + t18 * t15 - t20 * t17, t1 * r_i_i_C(1) + (t15 * t6 - t17 * t4) * r_i_i_C(2) + (-t14 * t15 + t17 * t21) * pkin(2), -t17 * t11, 0, 0, 0; 0, (-t19 * r_i_i_C(2) - t24) * t11, t13, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:18
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (112->37), mult. (271->57), div. (0->0), fcn. (353->10), ass. (0->28)
	t18 = sin(pkin(11));
	t21 = cos(pkin(11));
	t23 = sin(qJ(2));
	t25 = cos(qJ(2));
	t12 = t23 * t18 - t25 * t21;
	t37 = t25 * pkin(2);
	t22 = cos(pkin(6));
	t36 = t22 * t25;
	t19 = sin(pkin(6));
	t24 = sin(qJ(1));
	t34 = t24 * t19;
	t26 = cos(qJ(1));
	t32 = t26 * t19;
	t31 = r_i_i_C(3) + qJ(4);
	t29 = t25 * t18 + t23 * t21;
	t10 = t29 * t22;
	t3 = -t26 * t10 + t24 * t12;
	t30 = t24 * t10 + t26 * t12;
	t17 = sin(pkin(12));
	t20 = cos(pkin(12));
	t28 = t20 * r_i_i_C(1) - t17 * r_i_i_C(2) + pkin(3);
	t27 = t12 * t22;
	t16 = pkin(1) + t37;
	t11 = t22 * t23 * pkin(2) + (-pkin(8) - qJ(3)) * t19;
	t8 = t12 * t19;
	t5 = t24 * t27 - t26 * t29;
	t2 = -t24 * t29 - t26 * t27;
	t1 = [(t17 * t32 + t3 * t20) * r_i_i_C(1) + (-t3 * t17 + t20 * t32) * r_i_i_C(2) + t3 * pkin(3) - t24 * t16 - t26 * t11 + t31 * t2, -t31 * t30 + (-t23 * t26 - t24 * t36) * pkin(2) + t28 * t5, t34, -t5, 0, 0; (t17 * t34 - t20 * t30) * r_i_i_C(1) + (t17 * t30 + t20 * t34) * r_i_i_C(2) - t30 * pkin(3) + t26 * t16 - t24 * t11 - t31 * t5, -t31 * t3 + (-t23 * t24 + t26 * t36) * pkin(2) + t28 * t2, -t32, -t2, 0, 0; 0, -t28 * t8 + (t29 * t31 + t37) * t19, t22, t8, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (173->46), mult. (344->72), div. (0->0), fcn. (446->12), ass. (0->35)
	t25 = sin(pkin(11));
	t27 = cos(pkin(11));
	t30 = sin(qJ(2));
	t32 = cos(qJ(2));
	t14 = t30 * t25 - t32 * t27;
	t45 = pkin(4) * sin(pkin(12));
	t44 = t32 * pkin(2);
	t43 = r_i_i_C(3) + pkin(9) + qJ(4);
	t28 = cos(pkin(6));
	t42 = t28 * t32;
	t26 = sin(pkin(6));
	t31 = sin(qJ(1));
	t40 = t31 * t26;
	t33 = cos(qJ(1));
	t38 = t33 * t26;
	t36 = t32 * t25 + t30 * t27;
	t12 = t36 * t28;
	t5 = t33 * t12 - t31 * t14;
	t37 = t31 * t12 + t33 * t14;
	t19 = cos(pkin(12)) * pkin(4) + pkin(3);
	t23 = pkin(12) + qJ(5);
	t21 = sin(t23);
	t22 = cos(t23);
	t35 = t22 * r_i_i_C(1) - t21 * r_i_i_C(2) + t19;
	t34 = t14 * t28;
	t20 = pkin(1) + t44;
	t16 = t21 * t38;
	t13 = t28 * t30 * pkin(2) + (-pkin(8) - qJ(3)) * t26;
	t11 = t36 * t26;
	t10 = t14 * t26;
	t7 = t31 * t34 - t33 * t36;
	t4 = -t31 * t36 - t33 * t34;
	t2 = t21 * t40 - t22 * t37;
	t1 = t21 * t37 + t22 * t40;
	t3 = [t16 * r_i_i_C(1) - t31 * t20 - t35 * t5 + t43 * t4 + (-t13 + (t22 * r_i_i_C(2) + t45) * t26) * t33, -t43 * t37 + (-t30 * t33 - t31 * t42) * pkin(2) + t35 * t7, t40, -t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t37 * t19 + t33 * t20 - t43 * t7 + (t26 * t45 - t13) * t31, t43 * t5 + (-t30 * t31 + t33 * t42) * pkin(2) + t35 * t4, -t38, -t4, (-t5 * t21 - t22 * t38) * r_i_i_C(1) + (-t5 * t22 + t16) * r_i_i_C(2), 0; 0, -t35 * t10 + t43 * t11 + t26 * t44, t28, t10, (-t11 * t21 + t28 * t22) * r_i_i_C(1) + (-t11 * t22 - t28 * t21) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:39:18
	% EndTime: 2019-10-10 09:39:19
	% DurationCPUTime: 0.44s
	% Computational Cost: add. (384->64), mult. (758->96), div. (0->0), fcn. (998->14), ass. (0->45)
	t35 = sin(pkin(11));
	t40 = sin(qJ(2));
	t43 = cos(qJ(2));
	t54 = cos(pkin(11));
	t24 = -t43 * t35 - t40 * t54;
	t41 = sin(qJ(1));
	t44 = cos(qJ(1));
	t37 = cos(pkin(6));
	t47 = -t40 * t35 + t43 * t54;
	t46 = t47 * t37;
	t10 = t41 * t24 + t44 * t46;
	t39 = sin(qJ(6));
	t21 = t24 * t37;
	t11 = -t44 * t21 + t41 * t47;
	t33 = pkin(12) + qJ(5);
	t31 = sin(t33);
	t32 = cos(t33);
	t36 = sin(pkin(6));
	t55 = t44 * t36;
	t4 = t11 * t32 - t31 * t55;
	t42 = cos(qJ(6));
	t64 = t10 * t42 + t4 * t39;
	t63 = t10 * t39 - t4 * t42;
	t29 = cos(pkin(12)) * pkin(4) + pkin(3);
	t50 = t42 * r_i_i_C(1) - t39 * r_i_i_C(2) + pkin(5);
	t62 = pkin(10) + r_i_i_C(3);
	t45 = t62 * t31 + t50 * t32 + t29;
	t61 = t43 * pkin(2);
	t58 = t37 * t43;
	t56 = t41 * t36;
	t52 = -t37 * t40 * pkin(2) + (sin(pkin(12)) * pkin(4) + pkin(8) + qJ(3)) * t36;
	t51 = -t41 * t21 - t44 * t47;
	t38 = -pkin(9) - qJ(4);
	t49 = t39 * r_i_i_C(1) + t42 * r_i_i_C(2) - t38;
	t48 = -t11 * t31 - t32 * t55;
	t30 = pkin(1) + t61;
	t20 = t24 * t36;
	t19 = t47 * t36;
	t16 = -t20 * t32 + t37 * t31;
	t13 = t44 * t24 - t41 * t46;
	t8 = t31 * t56 - t32 * t51;
	t7 = -t31 * t51 - t32 * t56;
	t2 = -t13 * t39 + t8 * t42;
	t1 = -t13 * t42 - t8 * t39;
	t3 = [-t4 * pkin(5) + t63 * r_i_i_C(1) + t64 * r_i_i_C(2) - t10 * t38 - t11 * t29 - t41 * t30 + t52 * t44 + t62 * t48, (-t44 * t40 - t41 * t58) * pkin(2) - t49 * t51 + t45 * t13, t56, -t13, -t50 * t7 + t62 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t8 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t13 * t38 - t29 * t51 + t44 * t30 + t52 * t41 + t62 * t7, (-t41 * t40 + t44 * t58) * pkin(2) + t49 * t11 + t45 * t10, -t55, -t10, t62 * t4 + t50 * t48, -t64 * r_i_i_C(1) + t63 * r_i_i_C(2); 0, t45 * t19 - t49 * t20 + t36 * t61, t37, -t19, t62 * t16 + t50 * (t20 * t31 + t37 * t32), (-t16 * t39 - t19 * t42) * r_i_i_C(1) + (-t16 * t42 + t19 * t39) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end