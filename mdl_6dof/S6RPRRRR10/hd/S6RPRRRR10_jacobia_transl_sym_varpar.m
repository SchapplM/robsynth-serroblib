% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRR10
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR10_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR10_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR10_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
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
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (13->12), mult. (28->19), div. (0->0), fcn. (35->6), ass. (0->12)
	t1 = sin(pkin(13));
	t5 = sin(qJ(1));
	t11 = t5 * t1;
	t3 = cos(pkin(13));
	t10 = t5 * t3;
	t6 = cos(qJ(1));
	t9 = t6 * t1;
	t8 = t6 * t3;
	t2 = sin(pkin(6));
	t7 = t2 * (r_i_i_C(3) + qJ(2));
	t4 = cos(pkin(6));
	t12 = [(-t4 * t9 - t10) * r_i_i_C(1) + (-t4 * t8 + t11) * r_i_i_C(2) - t5 * pkin(1) + t6 * t7, t5 * t2, 0, 0, 0, 0; (-t11 * t4 + t8) * r_i_i_C(1) + (-t10 * t4 - t9) * r_i_i_C(2) + t6 * pkin(1) + t5 * t7, -t6 * t2, 0, 0, 0, 0; 0, t4, 0, 0, 0, 0;];
	Ja_transl = t12;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (62->32), mult. (166->58), div. (0->0), fcn. (213->10), ass. (0->29)
	t34 = r_i_i_C(3) + pkin(9);
	t13 = cos(pkin(7));
	t15 = sin(qJ(3));
	t17 = cos(qJ(3));
	t10 = sin(pkin(7));
	t11 = sin(pkin(6));
	t18 = cos(qJ(1));
	t26 = t18 * t11;
	t24 = t10 * t26;
	t12 = cos(pkin(13));
	t14 = cos(pkin(6));
	t29 = t14 * t18;
	t16 = sin(qJ(1));
	t9 = sin(pkin(13));
	t32 = t16 * t9;
	t3 = -t12 * t29 + t32;
	t27 = t16 * t12;
	t4 = t9 * t29 + t27;
	t33 = (t13 * t3 + t24) * t17 + t4 * t15;
	t30 = t13 * t15;
	t28 = t16 * t11;
	t25 = t11 * qJ(2);
	t5 = -t14 * t27 - t18 * t9;
	t22 = t10 * t28 + t13 * t5;
	t19 = t15 * t24 - t17 * t4 + t3 * t30;
	t6 = t12 * t18 - t14 * t32;
	t2 = t22 * t15 + t17 * t6;
	t1 = -t15 * t6 + t22 * t17;
	t7 = [t19 * r_i_i_C(1) + t33 * r_i_i_C(2) - t4 * pkin(2) - t16 * pkin(1) + t18 * t25 + t34 * (-t3 * t10 + t13 * t26), t28, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; t18 * pkin(1) + t6 * pkin(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t25 + t34 * (-t10 * t5 + t13 * t28), -t26, -t33 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0, 0; 0, t14, (t17 * r_i_i_C(1) - t15 * r_i_i_C(2)) * t14 * t10 + ((t12 * t13 * t17 - t15 * t9) * r_i_i_C(1) + (-t12 * t30 - t17 * t9) * r_i_i_C(2)) * t11, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (176->48), mult. (480->84), div. (0->0), fcn. (625->12), ass. (0->42)
	t27 = cos(pkin(6));
	t25 = cos(pkin(13));
	t33 = cos(qJ(1));
	t41 = t33 * t25;
	t22 = sin(pkin(13));
	t30 = sin(qJ(1));
	t46 = t30 * t22;
	t16 = -t27 * t41 + t46;
	t23 = sin(pkin(7));
	t26 = cos(pkin(7));
	t24 = sin(pkin(6));
	t42 = t33 * t24;
	t11 = -t16 * t23 + t26 * t42;
	t28 = sin(qJ(4));
	t31 = cos(qJ(4));
	t43 = t33 * t22;
	t44 = t30 * t25;
	t17 = t27 * t43 + t44;
	t29 = sin(qJ(3));
	t32 = cos(qJ(3));
	t39 = t23 * t42;
	t47 = t26 * t29;
	t6 = t16 * t47 - t17 * t32 + t29 * t39;
	t52 = t11 * t31 - t6 * t28;
	t51 = t11 * t28 + t6 * t31;
	t36 = t27 * t44 + t43;
	t45 = t30 * t24;
	t50 = -t23 * t45 + t36 * t26;
	t49 = r_i_i_C(3) + pkin(10);
	t48 = t23 * t27;
	t40 = t24 * qJ(2);
	t37 = t31 * r_i_i_C(1) - t28 * r_i_i_C(2) + pkin(3);
	t34 = -t17 * t29 + (-t16 * t26 - t39) * t32;
	t13 = t23 * t36 + t26 * t45;
	t18 = -t27 * t46 + t41;
	t15 = -t24 * t25 * t23 + t27 * t26;
	t10 = t29 * t48 + (t22 * t32 + t25 * t47) * t24;
	t8 = t18 * t32 - t50 * t29;
	t7 = t18 * t29 + t50 * t32;
	t2 = t13 * t28 + t8 * t31;
	t1 = t13 * t31 - t8 * t28;
	t3 = [-t30 * pkin(1) - t17 * pkin(2) + t6 * pkin(3) + t11 * pkin(9) + t51 * r_i_i_C(1) + t52 * r_i_i_C(2) + t33 * t40 + t49 * t34, t45, -t37 * t7 + t49 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t33 * pkin(1) + t18 * pkin(2) + t8 * pkin(3) + t13 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t30 * t40 + t49 * t7, -t42, t37 * t34 - t49 * t6, -t52 * r_i_i_C(1) + t51 * r_i_i_C(2), 0, 0; 0, t27, t49 * t10 + t37 * (t32 * t48 + (t25 * t26 * t32 - t22 * t29) * t24), (-t10 * t28 + t15 * t31) * r_i_i_C(1) + (-t10 * t31 - t15 * t28) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (290->57), mult. (669->95), div. (0->0), fcn. (869->14), ass. (0->48)
	t38 = sin(qJ(4));
	t66 = t38 * pkin(4) + pkin(9);
	t33 = sin(pkin(7));
	t36 = cos(pkin(7));
	t37 = cos(pkin(6));
	t32 = sin(pkin(13));
	t43 = cos(qJ(1));
	t54 = t43 * t32;
	t35 = cos(pkin(13));
	t40 = sin(qJ(1));
	t55 = t40 * t35;
	t47 = t37 * t55 + t54;
	t34 = sin(pkin(6));
	t56 = t40 * t34;
	t65 = -t33 * t56 + t47 * t36;
	t41 = cos(qJ(4));
	t28 = t41 * pkin(4) + pkin(3);
	t31 = qJ(4) + qJ(5);
	t29 = sin(t31);
	t30 = cos(t31);
	t48 = t30 * r_i_i_C(1) - t29 * r_i_i_C(2) + t28;
	t52 = t43 * t35;
	t57 = t40 * t32;
	t22 = -t37 * t52 + t57;
	t23 = t37 * t54 + t55;
	t39 = sin(qJ(3));
	t42 = cos(qJ(3));
	t53 = t43 * t34;
	t50 = t33 * t53;
	t58 = t36 * t39;
	t12 = t22 * t58 - t23 * t42 + t39 * t50;
	t17 = -t22 * t33 + t36 * t53;
	t64 = (t12 * t29 - t17 * t30) * r_i_i_C(1) + (t12 * t30 + t17 * t29) * r_i_i_C(2);
	t24 = -t37 * t57 + t52;
	t14 = t24 * t42 - t65 * t39;
	t19 = t47 * t33 + t36 * t56;
	t5 = -t14 * t29 + t19 * t30;
	t6 = t14 * t30 + t19 * t29;
	t63 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t59 = t33 * t37;
	t16 = t39 * t59 + (t32 * t42 + t35 * t58) * t34;
	t21 = -t34 * t35 * t33 + t37 * t36;
	t62 = (-t16 * t29 + t21 * t30) * r_i_i_C(1) + (-t16 * t30 - t21 * t29) * r_i_i_C(2);
	t60 = r_i_i_C(3) + pkin(11) + pkin(10);
	t51 = t34 * qJ(2);
	t45 = -t23 * t39 + (-t22 * t36 - t50) * t42;
	t13 = t24 * t39 + t65 * t42;
	t1 = [-t40 * pkin(1) - t23 * pkin(2) + t43 * t51 + t60 * t45 + t48 * t12 + (t29 * r_i_i_C(1) + t30 * r_i_i_C(2) + t66) * t17, t56, -t48 * t13 + t60 * t14, (-t14 * t38 + t19 * t41) * pkin(4) + t63, t63, 0; t43 * pkin(1) + t24 * pkin(2) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t60 * t13 + t14 * t28 + t66 * t19 + t40 * t51, -t53, -t12 * t60 + t48 * t45, (t12 * t38 - t17 * t41) * pkin(4) + t64, t64, 0; 0, t37, t60 * t16 + t48 * (t42 * t59 + (t35 * t36 * t42 - t32 * t39) * t34), (-t16 * t38 + t21 * t41) * pkin(4) + t62, t62, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (660->74), mult. (1502->121), div. (0->0), fcn. (1972->16), ass. (0->60)
	t52 = cos(qJ(1));
	t72 = sin(pkin(13));
	t76 = cos(pkin(6));
	t64 = t76 * t72;
	t74 = cos(pkin(13));
	t78 = sin(qJ(1));
	t35 = t52 * t64 + t78 * t74;
	t49 = sin(qJ(3));
	t79 = cos(qJ(3));
	t66 = t76 * t74;
	t34 = -t52 * t66 + t78 * t72;
	t46 = sin(pkin(6));
	t73 = sin(pkin(7));
	t68 = t46 * t73;
	t75 = cos(pkin(7));
	t87 = t34 * t75 + t52 * t68;
	t21 = t35 * t79 - t87 * t49;
	t69 = t46 * t75;
	t30 = t34 * t73 - t52 * t69;
	t45 = qJ(4) + qJ(5);
	t43 = sin(t45);
	t44 = cos(t45);
	t10 = t21 * t44 + t30 * t43;
	t20 = t35 * t49 + t87 * t79;
	t47 = sin(qJ(6));
	t50 = cos(qJ(6));
	t92 = t10 * t47 - t20 * t50;
	t91 = -t10 * t50 - t20 * t47;
	t48 = sin(qJ(4));
	t90 = pkin(4) * t48 + pkin(9);
	t83 = r_i_i_C(3) + pkin(12);
	t86 = t50 * r_i_i_C(1) - t47 * r_i_i_C(2) + pkin(5);
	t85 = -t21 * t43 + t30 * t44;
	t56 = t52 * t72 + t78 * t66;
	t84 = t56 * t75 - t78 * t68;
	t77 = t52 * t46;
	t71 = t78 * t46;
	t65 = t76 * t73;
	t63 = t75 * t74;
	t53 = -pkin(11) - pkin(10);
	t62 = t47 * r_i_i_C(1) + t50 * r_i_i_C(2) - t53;
	t60 = t83 * t10 + t86 * t85;
	t36 = t52 * t74 - t78 * t64;
	t25 = t36 * t79 - t84 * t49;
	t54 = t56 * t73 + t78 * t69;
	t13 = t25 * t43 - t54 * t44;
	t14 = t25 * t44 + t54 * t43;
	t59 = -t86 * t13 + t83 * t14;
	t28 = t49 * t65 + (t49 * t63 + t79 * t72) * t46;
	t33 = -t74 * t68 + t76 * t75;
	t19 = t28 * t44 + t33 * t43;
	t58 = t83 * t19 + t86 * (-t28 * t43 + t33 * t44);
	t51 = cos(qJ(4));
	t42 = t51 * pkin(4) + pkin(3);
	t57 = -t83 * t43 - t86 * t44 - t42;
	t27 = -t79 * t65 + (t49 * t72 - t63 * t79) * t46;
	t24 = t36 * t49 + t84 * t79;
	t2 = t14 * t50 + t24 * t47;
	t1 = -t14 * t47 + t24 * t50;
	t3 = [-t78 * pkin(1) - t35 * pkin(2) - t10 * pkin(5) + t91 * r_i_i_C(1) + t92 * r_i_i_C(2) + qJ(2) * t77 + t20 * t53 - t21 * t42 - t90 * t30 + t83 * t85, t71, t57 * t24 + t62 * t25, (-t25 * t48 + t54 * t51) * pkin(4) + t59, t59, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t52 * pkin(1) + t36 * pkin(2) + t14 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t71 + t83 * t13 - t24 * t53 + t25 * t42 + t90 * t54, -t77, t57 * t20 + t62 * t21, (-t21 * t48 + t30 * t51) * pkin(4) + t60, t60, -t92 * r_i_i_C(1) + t91 * r_i_i_C(2); 0, t76, t57 * t27 + t62 * t28, (-t28 * t48 + t33 * t51) * pkin(4) + t58, t58, (-t19 * t47 + t27 * t50) * r_i_i_C(1) + (-t19 * t50 - t27 * t47) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end