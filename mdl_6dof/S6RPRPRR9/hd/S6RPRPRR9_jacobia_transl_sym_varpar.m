% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR9_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR9_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR9_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:09
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
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:09
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (13->12), mult. (28->19), div. (0->0), fcn. (35->6), ass. (0->12)
	t1 = sin(pkin(12));
	t5 = sin(qJ(1));
	t11 = t5 * t1;
	t3 = cos(pkin(12));
	t10 = t5 * t3;
	t6 = cos(qJ(1));
	t9 = t6 * t1;
	t8 = t6 * t3;
	t2 = sin(pkin(6));
	t7 = t2 * (r_i_i_C(3) + qJ(2));
	t4 = cos(pkin(6));
	t12 = [(-t4 * t9 - t10) * r_i_i_C(1) + (-t4 * t8 + t11) * r_i_i_C(2) - t5 * pkin(1) + t6 * t7, t5 * t2, 0, 0, 0, 0; (-t4 * t11 + t8) * r_i_i_C(1) + (-t4 * t10 - t9) * r_i_i_C(2) + t6 * pkin(1) + t5 * t7, -t6 * t2, 0, 0, 0, 0; 0, t4, 0, 0, 0, 0;];
	Ja_transl = t12;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.18s
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
	t12 = cos(pkin(12));
	t14 = cos(pkin(6));
	t29 = t14 * t18;
	t16 = sin(qJ(1));
	t9 = sin(pkin(12));
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
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (120->56), mult. (315->91), div. (0->0), fcn. (404->12), ass. (0->38)
	t26 = sin(qJ(3));
	t41 = t26 * pkin(3);
	t28 = cos(qJ(3));
	t40 = t28 * pkin(3);
	t19 = sin(pkin(12));
	t27 = sin(qJ(1));
	t39 = t27 * t19;
	t21 = sin(pkin(6));
	t38 = t27 * t21;
	t23 = cos(pkin(12));
	t37 = t27 * t23;
	t29 = cos(qJ(1));
	t36 = t29 * t19;
	t35 = t29 * t21;
	t34 = t29 * t23;
	t33 = pkin(9) + qJ(4);
	t20 = sin(pkin(7));
	t24 = cos(pkin(7));
	t18 = sin(pkin(13));
	t22 = cos(pkin(13));
	t30 = t28 * t18 + t26 * t22;
	t4 = t30 * t20;
	t32 = r_i_i_C(1) * t4 + t20 * t41 + t33 * t24 + qJ(2);
	t25 = cos(pkin(6));
	t11 = -t25 * t37 - t36;
	t12 = -t25 * t39 + t34;
	t13 = t26 * t18 - t28 * t22;
	t6 = t30 * t24;
	t31 = t11 * t6 - t12 * t13;
	t17 = pkin(2) + t40;
	t10 = t25 * t36 + t37;
	t9 = -t25 * t34 + t39;
	t8 = -t33 * t20 + t24 * t41;
	t5 = t13 * t24;
	t3 = t13 * t20;
	t2 = -t11 * t20 + t24 * t38;
	t1 = -t11 * t5 - t12 * t30 - t3 * t38;
	t7 = [-t27 * pkin(1) + (t13 * r_i_i_C(1) + r_i_i_C(2) * t30 - t17) * t10 + (t6 * r_i_i_C(1) - t5 * r_i_i_C(2) - t20 * r_i_i_C(3) + t8) * t9 + (-r_i_i_C(2) * t3 + r_i_i_C(3) * t24 + t32) * t35, t38, t1 * r_i_i_C(1) + (-t4 * t38 - t31) * r_i_i_C(2) + (-t12 * t26 + (t11 * t24 + t20 * t38) * t28) * pkin(3), t2, 0, 0; t29 * pkin(1) + t31 * r_i_i_C(1) + t1 * r_i_i_C(2) + t2 * r_i_i_C(3) + t11 * t8 + t12 * t17 + t32 * t38, -t35, (-t10 * t30 + t3 * t35 + t9 * t5) * r_i_i_C(1) + (t10 * t13 + t4 * t35 + t9 * t6) * r_i_i_C(2) + (-t10 * t26 + (-t20 * t35 - t24 * t9) * t28) * pkin(3), t9 * t20 - t24 * t35, 0, 0; 0, t25, (-t3 * r_i_i_C(1) - t4 * r_i_i_C(2) + t20 * t40) * t25 + ((-t19 * t30 - t23 * t5) * r_i_i_C(1) + (t13 * t19 - t23 * t6) * r_i_i_C(2) + (t23 * t24 * t28 - t19 * t26) * pkin(3)) * t21, -t21 * t23 * t20 + t25 * t24, 0, 0;];
	Ja_transl = t7;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (297->66), mult. (797->111), div. (0->0), fcn. (1047->14), ass. (0->50)
	t40 = cos(pkin(6));
	t38 = cos(pkin(12));
	t46 = cos(qJ(1));
	t53 = t46 * t38;
	t34 = sin(pkin(12));
	t43 = sin(qJ(1));
	t58 = t43 * t34;
	t24 = -t40 * t53 + t58;
	t35 = sin(pkin(7));
	t39 = cos(pkin(7));
	t36 = sin(pkin(6));
	t54 = t46 * t36;
	t14 = -t24 * t35 + t39 * t54;
	t41 = sin(qJ(5));
	t44 = cos(qJ(5));
	t33 = sin(pkin(13));
	t37 = cos(pkin(13));
	t42 = sin(qJ(3));
	t45 = cos(qJ(3));
	t49 = t45 * t33 + t42 * t37;
	t18 = t49 * t35;
	t20 = t49 * t39;
	t55 = t46 * t34;
	t56 = t43 * t38;
	t25 = t40 * t55 + t56;
	t28 = t42 * t33 - t45 * t37;
	t7 = t18 * t54 + t24 * t20 + t25 * t28;
	t63 = t14 * t44 - t7 * t41;
	t62 = t14 * t41 + t7 * t44;
	t13 = t40 * t18 + (t20 * t38 - t28 * t34) * t36;
	t61 = -r_i_i_C(3) - pkin(10);
	t60 = pkin(3) * t42;
	t57 = t43 * t36;
	t52 = pkin(9) + qJ(4);
	t51 = t36 * (t35 * t60 + t52 * t39 + qJ(2));
	t48 = t44 * r_i_i_C(1) - t41 * r_i_i_C(2) + pkin(4);
	t17 = t28 * t35;
	t19 = t28 * t39;
	t47 = -t17 * t54 - t24 * t19 + t25 * t49;
	t26 = -t40 * t56 - t55;
	t27 = -t40 * t58 + t53;
	t10 = t18 * t57 + t26 * t20 - t27 * t28;
	t32 = t45 * pkin(3) + pkin(2);
	t23 = -t36 * t38 * t35 + t40 * t39;
	t22 = -t52 * t35 + t39 * t60;
	t16 = -t26 * t35 + t39 * t57;
	t9 = -t17 * t57 - t26 * t19 - t27 * t49;
	t2 = t10 * t44 + t16 * t41;
	t1 = -t10 * t41 + t16 * t44;
	t3 = [-t43 * pkin(1) + t7 * pkin(4) + t62 * r_i_i_C(1) + t63 * r_i_i_C(2) + t24 * t22 - t25 * t32 + t46 * t51 + t61 * t47, t57, -t61 * t10 + t48 * t9 + (-t27 * t42 + (t26 * t39 + t35 * t57) * t45) * pkin(3), t16, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t46 * pkin(1) + t10 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t26 * t22 + t27 * t32 + t43 * t51 + t61 * t9, -t54, t61 * t7 - t48 * t47 + (-t25 * t42 + (-t24 * t39 - t35 * t54) * t45) * pkin(3), -t14, -t63 * r_i_i_C(1) + t62 * r_i_i_C(2), 0; 0, t40, -t61 * t13 + t48 * (-t40 * t17 + (-t19 * t38 - t34 * t49) * t36) + (t35 * t40 * t45 + (t38 * t39 * t45 - t34 * t42) * t36) * pkin(3), t23, (-t13 * t41 + t23 * t44) * r_i_i_C(1) + (-t13 * t44 - t23 * t41) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (658->83), mult. (1787->137), div. (0->0), fcn. (2376->16), ass. (0->60)
	t42 = sin(pkin(13));
	t46 = cos(pkin(13));
	t52 = sin(qJ(3));
	t56 = cos(qJ(3));
	t37 = t52 * t42 - t56 * t46;
	t44 = sin(pkin(7));
	t26 = t37 * t44;
	t48 = cos(pkin(7));
	t28 = t37 * t48;
	t49 = cos(pkin(6));
	t47 = cos(pkin(12));
	t57 = cos(qJ(1));
	t66 = t57 * t47;
	t43 = sin(pkin(12));
	t53 = sin(qJ(1));
	t71 = t53 * t43;
	t33 = -t49 * t66 + t71;
	t68 = t57 * t43;
	t69 = t53 * t47;
	t34 = t49 * t68 + t69;
	t62 = t56 * t42 + t52 * t46;
	t45 = sin(pkin(6));
	t67 = t57 * t45;
	t14 = -t26 * t67 - t33 * t28 + t34 * t62;
	t50 = sin(qJ(6));
	t54 = cos(qJ(6));
	t27 = t62 * t44;
	t29 = t62 * t48;
	t15 = t27 * t67 + t33 * t29 + t34 * t37;
	t23 = -t33 * t44 + t48 * t67;
	t51 = sin(qJ(5));
	t55 = cos(qJ(5));
	t6 = t15 * t55 + t23 * t51;
	t79 = t14 * t54 + t6 * t50;
	t78 = -t14 * t50 + t6 * t54;
	t75 = t15 * t51 - t23 * t55;
	t21 = t49 * t27 + (t29 * t47 - t37 * t43) * t45;
	t74 = pkin(11) + r_i_i_C(3);
	t73 = pkin(3) * t52;
	t70 = t53 * t45;
	t65 = pkin(9) + qJ(4);
	t64 = t45 * (t44 * t73 + t65 * t48 + qJ(2));
	t61 = t54 * r_i_i_C(1) - t50 * r_i_i_C(2) + pkin(5);
	t60 = -t50 * r_i_i_C(1) - t54 * r_i_i_C(2) - pkin(10);
	t35 = -t49 * t69 - t68;
	t59 = -t35 * t44 + t48 * t70;
	t36 = -t49 * t71 + t66;
	t18 = t27 * t70 + t35 * t29 - t36 * t37;
	t58 = t74 * t51 + t61 * t55 + pkin(4);
	t41 = t56 * pkin(3) + pkin(2);
	t32 = -t45 * t47 * t44 + t49 * t48;
	t31 = -t65 * t44 + t48 * t73;
	t20 = -t49 * t26 + (-t28 * t47 - t43 * t62) * t45;
	t17 = -t26 * t70 - t35 * t28 - t36 * t62;
	t10 = t21 * t55 + t32 * t51;
	t8 = t18 * t55 + t59 * t51;
	t7 = t18 * t51 - t59 * t55;
	t2 = -t17 * t50 + t8 * t54;
	t1 = -t17 * t54 - t8 * t50;
	t3 = [-t53 * pkin(1) + t15 * pkin(4) + t6 * pkin(5) - t14 * pkin(10) + t78 * r_i_i_C(1) - t79 * r_i_i_C(2) + t33 * t31 - t34 * t41 + t57 * t64 + t74 * t75, t70, -t60 * t18 + (-t36 * t52 + (t35 * t48 + t44 * t70) * t56) * pkin(3) + t58 * t17, t59, -t61 * t7 + t74 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t57 * pkin(1) + t18 * pkin(4) + t8 * pkin(5) - t17 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t35 * t31 + t36 * t41 + t53 * t64 + t74 * t7, -t67, t60 * t15 + (-t34 * t52 + (-t33 * t48 - t44 * t67) * t56) * pkin(3) - t58 * t14, -t23, -t6 * t74 + t61 * t75, t79 * r_i_i_C(1) + t78 * r_i_i_C(2); 0, t49, -t60 * t21 + (t49 * t44 * t56 + (t47 * t48 * t56 - t43 * t52) * t45) * pkin(3) + t58 * t20, t32, t74 * t10 + t61 * (-t21 * t51 + t32 * t55), (-t10 * t50 - t20 * t54) * r_i_i_C(1) + (-t10 * t54 + t20 * t50) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end