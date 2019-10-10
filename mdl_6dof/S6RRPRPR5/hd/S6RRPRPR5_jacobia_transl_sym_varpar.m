% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:37
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:37
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:37
	% EndTime: 2019-10-10 10:11:38
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
	% StartTime: 2019-10-10 10:11:37
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.14s
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
	% StartTime: 2019-10-10 10:11:37
	% EndTime: 2019-10-10 10:11:38
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
	t19 = sin(pkin(11));
	t21 = cos(pkin(11));
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
	% StartTime: 2019-10-10 10:11:37
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (255->52), mult. (641->83), div. (0->0), fcn. (844->12), ass. (0->37)
	t38 = cos(qJ(2));
	t51 = t38 * pkin(2);
	t33 = cos(pkin(6));
	t50 = t33 * t38;
	t30 = sin(pkin(6));
	t36 = sin(qJ(1));
	t49 = t36 * t30;
	t39 = cos(qJ(1));
	t48 = t39 * t30;
	t47 = r_i_i_C(3) + qJ(5);
	t34 = sin(qJ(4));
	t37 = cos(qJ(4));
	t29 = sin(pkin(11));
	t32 = cos(pkin(11));
	t35 = sin(qJ(2));
	t44 = t38 * t29 + t35 * t32;
	t19 = t44 * t33;
	t21 = t35 * t29 - t38 * t32;
	t9 = t39 * t19 - t36 * t21;
	t46 = -t34 * t48 + t9 * t37;
	t45 = t36 * t19 + t39 * t21;
	t28 = sin(pkin(12));
	t31 = cos(pkin(12));
	t43 = r_i_i_C(1) * t31 - r_i_i_C(2) * t28 + pkin(4);
	t42 = t28 * r_i_i_C(1) + t31 * r_i_i_C(2) + pkin(9);
	t1 = t9 * t34 + t37 * t48;
	t41 = t21 * t33;
	t40 = t47 * t34 + t43 * t37 + pkin(3);
	t27 = pkin(1) + t51;
	t20 = t33 * t35 * pkin(2) + (-pkin(8) - qJ(3)) * t30;
	t18 = t44 * t30;
	t13 = t18 * t34 - t33 * t37;
	t11 = t36 * t41 - t39 * t44;
	t8 = -t36 * t44 - t39 * t41;
	t6 = t34 * t49 - t37 * t45;
	t5 = -t34 * t45 - t37 * t49;
	t2 = [(t8 * t28 - t31 * t46) * r_i_i_C(1) + (t28 * t46 + t8 * t31) * r_i_i_C(2) - t46 * pkin(4) - t9 * pkin(3) + t8 * pkin(9) - t36 * t27 - t39 * t20 - t47 * t1, (-t39 * t35 - t36 * t50) * pkin(2) - t42 * t45 + t40 * t11, t49, -t43 * t5 + t47 * t6, t5, 0; (-t11 * t28 + t6 * t31) * r_i_i_C(1) + (-t11 * t31 - t6 * t28) * r_i_i_C(2) + t6 * pkin(4) - t45 * pkin(3) - t11 * pkin(9) + t39 * t27 - t36 * t20 + t47 * t5, (-t36 * t35 + t39 * t50) * pkin(2) + t42 * t9 + t40 * t8, -t48, -t43 * t1 + t47 * t46, t1, 0; 0, t42 * t18 + (-t21 * t40 + t51) * t30, t33, t47 * (t18 * t37 + t33 * t34) - t43 * t13, t13, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:11:37
	% EndTime: 2019-10-10 10:11:38
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (370->58), mult. (809->94), div. (0->0), fcn. (1065->14), ass. (0->44)
	t36 = sin(pkin(11));
	t41 = sin(qJ(2));
	t44 = cos(qJ(2));
	t54 = cos(pkin(11));
	t49 = -t41 * t36 + t44 * t54;
	t40 = sin(qJ(4));
	t43 = cos(qJ(4));
	t30 = cos(pkin(12)) * pkin(5) + pkin(4);
	t34 = pkin(12) + qJ(6);
	t32 = sin(t34);
	t33 = cos(t34);
	t50 = t33 * r_i_i_C(1) - t32 * r_i_i_C(2) + t30;
	t59 = r_i_i_C(3) + pkin(10) + qJ(5);
	t46 = t59 * t40 + t50 * t43 + pkin(3);
	t60 = t44 * pkin(2);
	t38 = cos(pkin(6));
	t58 = t38 * t44;
	t37 = sin(pkin(6));
	t42 = sin(qJ(1));
	t56 = t42 * t37;
	t45 = cos(qJ(1));
	t55 = t45 * t37;
	t53 = -sin(pkin(12)) * pkin(5) - pkin(9);
	t24 = -t44 * t36 - t41 * t54;
	t21 = t24 * t38;
	t11 = -t45 * t21 + t42 * t49;
	t4 = t11 * t43 - t40 * t55;
	t51 = -t42 * t21 - t45 * t49;
	t3 = t11 * t40 + t43 * t55;
	t48 = t32 * r_i_i_C(1) + t33 * r_i_i_C(2) - t53;
	t47 = t49 * t38;
	t31 = pkin(1) + t60;
	t22 = t38 * t41 * pkin(2) + (-pkin(8) - qJ(3)) * t37;
	t20 = t24 * t37;
	t19 = t49 * t37;
	t16 = -t20 * t43 + t38 * t40;
	t15 = -t20 * t40 - t38 * t43;
	t13 = t45 * t24 - t42 * t47;
	t10 = t42 * t24 + t45 * t47;
	t8 = t40 * t56 - t43 * t51;
	t7 = -t40 * t51 - t43 * t56;
	t2 = -t13 * t32 + t8 * t33;
	t1 = -t13 * t33 - t8 * t32;
	t5 = [-t11 * pkin(3) + t48 * t10 - t45 * t22 - t59 * t3 - t42 * t31 - t50 * t4, (-t45 * t41 - t42 * t58) * pkin(2) - t48 * t51 + t46 * t13, t56, -t50 * t7 + t59 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -pkin(3) * t51 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t53 * t13 - t42 * t22 + t8 * t30 + t45 * t31 + t59 * t7, (-t42 * t41 + t45 * t58) * pkin(2) + t48 * t11 + t46 * t10, -t55, -t50 * t3 + t59 * t4, t3, (-t10 * t33 - t4 * t32) * r_i_i_C(1) + (t10 * t32 - t4 * t33) * r_i_i_C(2); 0, t46 * t19 - t48 * t20 + t37 * t60, t38, -t50 * t15 + t59 * t16, t15, (-t16 * t32 - t19 * t33) * r_i_i_C(1) + (-t16 * t33 + t19 * t32) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end