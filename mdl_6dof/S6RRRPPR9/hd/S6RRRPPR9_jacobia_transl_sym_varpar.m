% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR9_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR9_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR9_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
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
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
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
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.15s
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
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (132->41), mult. (334->69), div. (0->0), fcn. (423->10), ass. (0->29)
	t20 = sin(qJ(3));
	t23 = cos(qJ(3));
	t17 = sin(pkin(11));
	t19 = cos(pkin(11));
	t28 = r_i_i_C(1) * t19 - r_i_i_C(2) * t17 + pkin(3);
	t33 = r_i_i_C(3) + qJ(4);
	t37 = t33 * t20 + t28 * t23 + pkin(2);
	t18 = sin(pkin(6));
	t22 = sin(qJ(1));
	t36 = t18 * t22;
	t35 = t18 * t23;
	t25 = cos(qJ(1));
	t34 = t18 * t25;
	t32 = cos(pkin(6));
	t21 = sin(qJ(2));
	t24 = cos(qJ(2));
	t29 = t25 * t32;
	t10 = t21 * t29 + t22 * t24;
	t31 = t10 * t23 - t20 * t34;
	t30 = t22 * t32;
	t27 = t17 * r_i_i_C(1) + t19 * r_i_i_C(2) + pkin(9);
	t1 = t10 * t20 + t23 * t34;
	t12 = -t21 * t30 + t25 * t24;
	t11 = t25 * t21 + t24 * t30;
	t9 = t22 * t21 - t24 * t29;
	t7 = t18 * t21 * t20 - t32 * t23;
	t6 = t12 * t23 + t20 * t36;
	t5 = t12 * t20 - t22 * t35;
	t2 = [(-t9 * t17 - t19 * t31) * r_i_i_C(1) + (t17 * t31 - t9 * t19) * r_i_i_C(2) - t31 * pkin(3) - t10 * pkin(2) - t9 * pkin(9) - t22 * pkin(1) + pkin(8) * t34 - t33 * t1, -t11 * t37 + t27 * t12, -t28 * t5 + t33 * t6, t5, 0, 0; (t11 * t17 + t6 * t19) * r_i_i_C(1) + (t11 * t19 - t6 * t17) * r_i_i_C(2) + t6 * pkin(3) + t12 * pkin(2) + t11 * pkin(9) + t25 * pkin(1) + pkin(8) * t36 + t33 * t5, t27 * t10 - t37 * t9, -t28 * t1 + t33 * t31, t1, 0, 0; 0, (t27 * t21 + t37 * t24) * t18, t33 * (t32 * t20 + t21 * t35) - t28 * t7, t7, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (198->54), mult. (504->91), div. (0->0), fcn. (646->10), ass. (0->36)
	t31 = sin(qJ(2));
	t32 = sin(qJ(1));
	t34 = cos(qJ(2));
	t35 = cos(qJ(1));
	t41 = cos(pkin(6));
	t38 = t35 * t41;
	t20 = t31 * t38 + t32 * t34;
	t30 = sin(qJ(3));
	t33 = cos(qJ(3));
	t28 = sin(pkin(6));
	t46 = t28 * t35;
	t10 = t20 * t33 - t30 * t46;
	t19 = t32 * t31 - t34 * t38;
	t27 = sin(pkin(11));
	t29 = cos(pkin(11));
	t53 = t10 * t27 - t19 * t29;
	t43 = r_i_i_C(2) + qJ(4);
	t52 = pkin(3) * t33 + t43 * t30 + pkin(2);
	t51 = pkin(4) + r_i_i_C(1);
	t49 = t27 * t33;
	t48 = t28 * t32;
	t47 = t28 * t33;
	t45 = t29 * t33;
	t44 = t33 * t34;
	t42 = r_i_i_C(3) + qJ(5);
	t39 = t32 * t41;
	t9 = t20 * t30 + t33 * t46;
	t36 = -t42 * t27 - t51 * t29 - pkin(3);
	t22 = -t31 * t39 + t35 * t34;
	t21 = t35 * t31 + t34 * t39;
	t18 = t41 * t30 + t31 * t47;
	t17 = t28 * t31 * t30 - t41 * t33;
	t14 = t22 * t33 + t30 * t48;
	t13 = t22 * t30 - t32 * t47;
	t3 = t14 * t27 - t21 * t29;
	t1 = [pkin(8) * t46 - t32 * pkin(1) - t20 * pkin(2) - t10 * pkin(3) - t19 * pkin(9) + t51 * (-t10 * t29 - t19 * t27) - t43 * t9 - t42 * t53, t22 * pkin(9) + t51 * (-t21 * t45 + t22 * t27) + t42 * (-t21 * t49 - t22 * t29) - t52 * t21, t36 * t13 + t43 * t14, t13, t3, 0; pkin(8) * t48 + t35 * pkin(1) + t22 * pkin(2) + t14 * pkin(3) + t21 * pkin(9) + t51 * (t14 * t29 + t21 * t27) + t42 * t3 + t43 * t13, t20 * pkin(9) + t51 * (-t19 * t45 + t20 * t27) + t42 * (-t19 * t49 - t20 * t29) - t52 * t19, t43 * t10 + t36 * t9, t9, t53, 0; 0, (t51 * (t27 * t31 + t29 * t44) + t42 * (t27 * t44 - t29 * t31) + pkin(9) * t31 + t52 * t34) * t28, t36 * t17 + t43 * t18, t17, t28 * t34 * t29 + t18 * t27, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (350->69), mult. (908->117), div. (0->0), fcn. (1183->12), ass. (0->46)
	t38 = sin(qJ(2));
	t40 = cos(qJ(2));
	t55 = cos(pkin(6));
	t61 = cos(qJ(1));
	t46 = t55 * t61;
	t59 = sin(qJ(1));
	t26 = t38 * t46 + t59 * t40;
	t37 = sin(qJ(3));
	t34 = sin(pkin(6));
	t52 = t34 * t61;
	t60 = cos(qJ(3));
	t16 = t26 * t60 - t37 * t52;
	t25 = t59 * t38 - t40 * t46;
	t33 = sin(pkin(11));
	t35 = cos(pkin(11));
	t3 = t16 * t33 - t25 * t35;
	t4 = t16 * t35 + t25 * t33;
	t54 = -r_i_i_C(3) - pkin(10) + qJ(4);
	t63 = t60 * pkin(3) + t54 * t37 + pkin(2);
	t62 = pkin(4) + pkin(5);
	t56 = t34 * t40;
	t53 = t33 * t60;
	t51 = t34 * t60;
	t50 = t34 * t59;
	t49 = t35 * t60;
	t48 = t40 * t60;
	t45 = t55 * t59;
	t36 = sin(qJ(6));
	t39 = cos(qJ(6));
	t44 = t36 * r_i_i_C(1) + t39 * r_i_i_C(2) + qJ(5);
	t43 = t39 * r_i_i_C(1) - t36 * r_i_i_C(2) + t62;
	t15 = t26 * t37 + t61 * t51;
	t41 = -t44 * t33 - t43 * t35 - pkin(3);
	t28 = -t38 * t45 + t61 * t40;
	t27 = t61 * t38 + t40 * t45;
	t24 = t55 * t37 + t38 * t51;
	t23 = t34 * t38 * t37 - t55 * t60;
	t20 = t28 * t60 + t37 * t50;
	t19 = t28 * t37 - t60 * t50;
	t14 = t24 * t35 - t33 * t56;
	t13 = t24 * t33 + t35 * t56;
	t8 = t20 * t35 + t27 * t33;
	t7 = t20 * t33 - t27 * t35;
	t2 = t7 * t36 + t8 * t39;
	t1 = -t8 * t36 + t7 * t39;
	t5 = [-t59 * pkin(1) - t26 * pkin(2) - t16 * pkin(3) + pkin(8) * t52 - t25 * pkin(9) - t54 * t15 - t44 * t3 - t43 * t4, t28 * pkin(9) + t44 * (-t27 * t53 - t28 * t35) + t43 * (-t27 * t49 + t28 * t33) - t63 * t27, t41 * t19 + t54 * t20, t19, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t61 * pkin(1) + t28 * pkin(2) + t20 * pkin(3) + pkin(8) * t50 + t27 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(5) + t54 * t19 + t62 * t8, t26 * pkin(9) + t44 * (-t25 * t53 - t26 * t35) + t43 * (-t25 * t49 + t26 * t33) - t63 * t25, t41 * t15 + t54 * t16, t15, t3, (t3 * t39 - t4 * t36) * r_i_i_C(1) + (-t3 * t36 - t4 * t39) * r_i_i_C(2); 0, (t44 * (t33 * t48 - t35 * t38) + t43 * (t33 * t38 + t35 * t48) + t38 * pkin(9) + t63 * t40) * t34, t41 * t23 + t54 * t24, t23, t13, (t13 * t39 - t14 * t36) * r_i_i_C(1) + (-t13 * t36 - t14 * t39) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end