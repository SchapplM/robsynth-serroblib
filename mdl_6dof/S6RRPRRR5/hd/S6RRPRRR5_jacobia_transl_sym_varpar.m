% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:57
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:44
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:44
	% EndTime: 2019-10-10 10:57:44
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
	% StartTime: 2019-10-10 10:57:44
	% EndTime: 2019-10-10 10:57:45
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
	% StartTime: 2019-10-10 10:57:44
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (49->24), mult. (114->41), div. (0->0), fcn. (143->8), ass. (0->19)
	t16 = cos(qJ(2));
	t10 = sin(pkin(12));
	t12 = cos(pkin(12));
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
	% StartTime: 2019-10-10 10:57:44
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (128->40), mult. (313->67), div. (0->0), fcn. (405->10), ass. (0->32)
	t38 = r_i_i_C(3) + pkin(9);
	t27 = cos(qJ(2));
	t37 = t27 * pkin(2);
	t22 = cos(pkin(6));
	t36 = t22 * t27;
	t20 = sin(pkin(6));
	t25 = sin(qJ(1));
	t35 = t25 * t20;
	t28 = cos(qJ(1));
	t34 = t28 * t20;
	t23 = sin(qJ(4));
	t26 = cos(qJ(4));
	t19 = sin(pkin(12));
	t21 = cos(pkin(12));
	t24 = sin(qJ(2));
	t31 = t27 * t19 + t24 * t21;
	t12 = t31 * t22;
	t14 = t24 * t19 - t27 * t21;
	t5 = t28 * t12 - t25 * t14;
	t33 = t23 * t34 - t5 * t26;
	t32 = t25 * t12 + t28 * t14;
	t30 = t26 * r_i_i_C(1) - t23 * r_i_i_C(2) + pkin(3);
	t29 = t5 * t23 + t26 * t34;
	t18 = pkin(1) + t37;
	t13 = t22 * t24 * pkin(2) + (-pkin(8) - qJ(3)) * t20;
	t11 = t14 * t22;
	t10 = t31 * t20;
	t7 = t25 * t11 - t28 * t31;
	t4 = -t28 * t11 - t25 * t31;
	t2 = t23 * t35 - t26 * t32;
	t1 = t23 * t32 + t26 * t35;
	t3 = [-t5 * pkin(3) + t33 * r_i_i_C(1) + t29 * r_i_i_C(2) - t28 * t13 - t25 * t18 + t38 * t4, -t38 * t32 + (-t24 * t28 - t25 * t36) * pkin(2) + t30 * t7, t35, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; -pkin(3) * t32 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t25 * t13 + t28 * t18 - t38 * t7, t38 * t5 + (-t24 * t25 + t28 * t36) * pkin(2) + t30 * t4, -t34, -t29 * r_i_i_C(1) + t33 * r_i_i_C(2), 0, 0; 0, t38 * t10 + (-t14 * t30 + t37) * t20, t22, (-t10 * t23 + t22 * t26) * r_i_i_C(1) + (-t10 * t26 - t22 * t23) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:44
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (288->57), mult. (727->94), div. (0->0), fcn. (957->12), ass. (0->42)
	t30 = sin(pkin(12));
	t35 = sin(qJ(2));
	t39 = cos(qJ(2));
	t49 = cos(pkin(12));
	t24 = -t30 * t39 - t35 * t49;
	t36 = sin(qJ(1));
	t40 = cos(qJ(1));
	t32 = cos(pkin(6));
	t43 = -t35 * t30 + t39 * t49;
	t42 = t43 * t32;
	t10 = t36 * t24 + t40 * t42;
	t33 = sin(qJ(5));
	t37 = cos(qJ(5));
	t21 = t24 * t32;
	t11 = -t21 * t40 + t36 * t43;
	t34 = sin(qJ(4));
	t38 = cos(qJ(4));
	t31 = sin(pkin(6));
	t50 = t40 * t31;
	t4 = t11 * t38 - t34 * t50;
	t59 = t10 * t37 + t33 * t4;
	t58 = t10 * t33 - t37 * t4;
	t46 = r_i_i_C(1) * t37 - r_i_i_C(2) * t33 + pkin(4);
	t57 = r_i_i_C(3) + pkin(10);
	t41 = t34 * t57 + t38 * t46 + pkin(3);
	t56 = t39 * pkin(2);
	t53 = t32 * t39;
	t51 = t36 * t31;
	t47 = -t21 * t36 - t40 * t43;
	t45 = r_i_i_C(1) * t33 + r_i_i_C(2) * t37 + pkin(9);
	t44 = -t11 * t34 - t38 * t50;
	t29 = pkin(1) + t56;
	t22 = t32 * t35 * pkin(2) + (-pkin(8) - qJ(3)) * t31;
	t20 = t24 * t31;
	t19 = t43 * t31;
	t16 = -t20 * t38 + t32 * t34;
	t13 = t24 * t40 - t36 * t42;
	t8 = t34 * t51 - t38 * t47;
	t7 = -t34 * t47 - t38 * t51;
	t2 = -t13 * t33 + t37 * t8;
	t1 = -t13 * t37 - t33 * t8;
	t3 = [-t11 * pkin(3) - t4 * pkin(4) + t10 * pkin(9) + r_i_i_C(1) * t58 + r_i_i_C(2) * t59 - t40 * t22 - t36 * t29 + t57 * t44, (-t35 * t40 - t36 * t53) * pkin(2) - t45 * t47 + t41 * t13, t51, -t46 * t7 + t57 * t8, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; -pkin(3) * t47 + t8 * pkin(4) - t13 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t36 * t22 + t40 * t29 + t57 * t7, (-t35 * t36 + t40 * t53) * pkin(2) + t45 * t11 + t41 * t10, -t50, t4 * t57 + t44 * t46, -r_i_i_C(1) * t59 + r_i_i_C(2) * t58, 0; 0, t19 * t41 - t45 * t20 + t31 * t56, t32, t57 * t16 + t46 * (t20 * t34 + t32 * t38), (-t16 * t33 - t19 * t37) * r_i_i_C(1) + (-t16 * t37 + t19 * t33) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:57:44
	% EndTime: 2019-10-10 10:57:45
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (438->64), mult. (955->103), div. (0->0), fcn. (1255->14), ass. (0->48)
	t40 = sin(pkin(12));
	t45 = sin(qJ(2));
	t49 = cos(qJ(2));
	t61 = cos(pkin(12));
	t55 = -t45 * t40 + t49 * t61;
	t44 = sin(qJ(4));
	t48 = cos(qJ(4));
	t47 = cos(qJ(5));
	t35 = t47 * pkin(5) + pkin(4);
	t39 = qJ(5) + qJ(6);
	t37 = sin(t39);
	t38 = cos(t39);
	t57 = r_i_i_C(1) * t38 - r_i_i_C(2) * t37 + t35;
	t66 = r_i_i_C(3) + pkin(11) + pkin(10);
	t52 = t66 * t44 + t57 * t48 + pkin(3);
	t30 = -t49 * t40 - t45 * t61;
	t42 = cos(pkin(6));
	t27 = t30 * t42;
	t46 = sin(qJ(1));
	t50 = cos(qJ(1));
	t17 = -t50 * t27 + t46 * t55;
	t41 = sin(pkin(6));
	t62 = t50 * t41;
	t10 = t17 * t48 - t44 * t62;
	t53 = t55 * t42;
	t16 = t46 * t30 + t50 * t53;
	t70 = (-t10 * t37 - t16 * t38) * r_i_i_C(1) + (-t10 * t38 + t16 * t37) * r_i_i_C(2);
	t58 = -t46 * t27 - t50 * t55;
	t63 = t46 * t41;
	t14 = t44 * t63 - t48 * t58;
	t19 = t50 * t30 - t46 * t53;
	t5 = -t14 * t37 - t19 * t38;
	t6 = t14 * t38 - t19 * t37;
	t69 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t26 = t30 * t41;
	t22 = -t26 * t48 + t42 * t44;
	t25 = t55 * t41;
	t68 = (-t22 * t37 - t25 * t38) * r_i_i_C(1) + (-t22 * t38 + t25 * t37) * r_i_i_C(2);
	t67 = t49 * pkin(2);
	t65 = t42 * t49;
	t43 = sin(qJ(5));
	t60 = -t43 * pkin(5) - pkin(9);
	t56 = -t17 * t44 - t48 * t62;
	t54 = t37 * r_i_i_C(1) + t38 * r_i_i_C(2) - t60;
	t36 = pkin(1) + t67;
	t28 = t42 * t45 * pkin(2) + (-pkin(8) - qJ(3)) * t41;
	t13 = -t44 * t58 - t48 * t63;
	t1 = [-t17 * pkin(3) - t57 * t10 + t54 * t16 - t50 * t28 - t46 * t36 + t66 * t56, (-t50 * t45 - t46 * t65) * pkin(2) - t54 * t58 + t52 * t19, t63, -t57 * t13 + t66 * t14, (-t14 * t43 - t19 * t47) * pkin(5) + t69, t69; -pkin(3) * t58 + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t66 * t13 + t14 * t35 + t60 * t19 - t46 * t28 + t50 * t36, (-t46 * t45 + t50 * t65) * pkin(2) + t54 * t17 + t52 * t16, -t62, t66 * t10 + t57 * t56, (-t10 * t43 - t16 * t47) * pkin(5) + t70, t70; 0, t52 * t25 - t54 * t26 + t41 * t67, t42, t66 * t22 + t57 * (t26 * t44 + t42 * t48), (-t22 * t43 - t25 * t47) * pkin(5) + t68, t68;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end