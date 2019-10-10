% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR9_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR9_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR9_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobia_transl_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
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
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.10s
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
	t12 = [(-t4 * t9 - t10) * r_i_i_C(1) + (-t4 * t8 + t11) * r_i_i_C(2) - t5 * pkin(1) + t6 * t7, t5 * t2, 0, 0, 0, 0; (-t11 * t4 + t8) * r_i_i_C(1) + (-t10 * t4 - t9) * r_i_i_C(2) + t6 * pkin(1) + t5 * t7, -t6 * t2, 0, 0, 0, 0; 0, t4, 0, 0, 0, 0;];
	Ja_transl = t12;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
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
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (176->48), mult. (480->84), div. (0->0), fcn. (625->12), ass. (0->42)
	t27 = cos(pkin(6));
	t25 = cos(pkin(12));
	t33 = cos(qJ(1));
	t41 = t33 * t25;
	t22 = sin(pkin(12));
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
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (249->58), mult. (590->97), div. (0->0), fcn. (767->14), ass. (0->49)
	t35 = sin(qJ(4));
	t60 = t35 * pkin(4) + pkin(9);
	t38 = cos(qJ(4));
	t24 = t38 * pkin(4) + pkin(3);
	t27 = qJ(4) + pkin(13);
	t25 = sin(t27);
	t26 = cos(t27);
	t59 = t26 * r_i_i_C(1) - t25 * r_i_i_C(2) + t24;
	t57 = r_i_i_C(3) + qJ(5) + pkin(10);
	t56 = cos(qJ(3));
	t29 = sin(pkin(7));
	t36 = sin(qJ(3));
	t55 = t29 * t36;
	t30 = sin(pkin(6));
	t31 = cos(pkin(12));
	t54 = t30 * t31;
	t32 = cos(pkin(7));
	t53 = t32 * t36;
	t28 = sin(pkin(12));
	t37 = sin(qJ(1));
	t52 = t37 * t28;
	t51 = t37 * t30;
	t50 = t37 * t31;
	t39 = cos(qJ(1));
	t49 = t39 * t28;
	t48 = t39 * t30;
	t47 = t39 * t31;
	t46 = t30 * qJ(2);
	t45 = t29 * t56;
	t44 = t32 * t56;
	t43 = t30 * t45;
	t33 = cos(pkin(6));
	t17 = -t33 * t47 + t52;
	t11 = -t17 * t29 + t32 * t48;
	t41 = t33 * t50 + t49;
	t40 = t41 * t32;
	t18 = t33 * t49 + t50;
	t4 = -t17 * t53 + t18 * t56 - t48 * t55;
	t13 = t41 * t29 + t32 * t51;
	t3 = t17 * t44 + t18 * t36 + t39 * t43;
	t19 = -t33 * t52 + t47;
	t16 = -t29 * t54 + t33 * t32;
	t10 = t33 * t55 + (t56 * t28 + t31 * t53) * t30;
	t9 = t30 * t28 * t36 - t33 * t45 - t44 * t54;
	t8 = t19 * t56 + (t29 * t51 - t40) * t36;
	t7 = t19 * t36 - t37 * t43 + t56 * t40;
	t2 = t13 * t25 + t8 * t26;
	t1 = t13 * t26 - t8 * t25;
	t5 = [-t37 * pkin(1) - t18 * pkin(2) - t57 * t3 + t39 * t46 - t59 * t4 + (t25 * r_i_i_C(1) + t26 * r_i_i_C(2) + t60) * t11, t51, t57 * t8 - t59 * t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + (t13 * t38 - t35 * t8) * pkin(4), t7, 0; t39 * pkin(1) + t19 * pkin(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t60 * t13 + t8 * t24 + t37 * t46 + t57 * t7, -t48, -t3 * t59 + t57 * t4, (-t11 * t26 - t4 * t25) * r_i_i_C(1) + (t11 * t25 - t4 * t26) * r_i_i_C(2) + (-t11 * t38 - t4 * t35) * pkin(4), t3, 0; 0, t33, t57 * t10 - t59 * t9, (-t10 * t25 + t16 * t26) * r_i_i_C(1) + (-t10 * t26 - t16 * t25) * r_i_i_C(2) + (-t10 * t35 + t16 * t38) * pkin(4), t9, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (538->75), mult. (1252->124), div. (0->0), fcn. (1645->16), ass. (0->57)
	t45 = cos(qJ(1));
	t61 = sin(pkin(12));
	t64 = cos(pkin(6));
	t53 = t64 * t61;
	t62 = cos(pkin(12));
	t66 = sin(qJ(1));
	t26 = t45 * t53 + t66 * t62;
	t37 = sin(pkin(7));
	t42 = sin(qJ(3));
	t54 = t64 * t62;
	t25 = -t45 * t54 + t66 * t61;
	t63 = cos(pkin(7));
	t58 = t25 * t63;
	t38 = sin(pkin(6));
	t67 = cos(qJ(3));
	t60 = t38 * t67;
	t11 = t45 * t37 * t60 + t26 * t42 + t67 * t58;
	t65 = t45 * t38;
	t12 = t26 * t67 + (-t37 * t65 - t58) * t42;
	t57 = t38 * t63;
	t21 = t25 * t37 - t45 * t57;
	t36 = qJ(4) + pkin(13);
	t34 = sin(t36);
	t35 = cos(t36);
	t4 = t12 * t35 + t21 * t34;
	t40 = sin(qJ(6));
	t43 = cos(qJ(6));
	t76 = -t11 * t43 + t4 * t40;
	t75 = -t11 * t40 - t4 * t43;
	t41 = sin(qJ(4));
	t74 = pkin(4) * t41 + pkin(9);
	t71 = -t12 * t34 + t21 * t35;
	t48 = t45 * t61 + t66 * t54;
	t59 = t66 * t38;
	t70 = -t37 * t59 + t48 * t63;
	t69 = r_i_i_C(3) + pkin(11);
	t56 = t64 * t37;
	t52 = t63 * t62;
	t51 = t43 * r_i_i_C(1) - t40 * r_i_i_C(2) + pkin(5);
	t39 = -qJ(5) - pkin(10);
	t50 = t40 * r_i_i_C(1) + t43 * r_i_i_C(2) - t39;
	t44 = cos(qJ(4));
	t33 = t44 * pkin(4) + pkin(3);
	t49 = -t69 * t34 - t51 * t35 - t33;
	t46 = t48 * t37 + t66 * t57;
	t27 = t45 * t62 - t66 * t53;
	t24 = -t38 * t62 * t37 + t64 * t63;
	t19 = t42 * t56 + (t42 * t52 + t67 * t61) * t38;
	t18 = t38 * t61 * t42 - t52 * t60 - t67 * t56;
	t16 = t27 * t67 - t70 * t42;
	t15 = t27 * t42 + t70 * t67;
	t10 = t19 * t35 + t24 * t34;
	t8 = t16 * t35 + t34 * t46;
	t7 = t16 * t34 - t35 * t46;
	t2 = t15 * t40 + t8 * t43;
	t1 = t15 * t43 - t8 * t40;
	t3 = [-t66 * pkin(1) - t26 * pkin(2) - t4 * pkin(5) + t75 * r_i_i_C(1) + t76 * r_i_i_C(2) + qJ(2) * t65 + t11 * t39 - t12 * t33 - t74 * t21 + t69 * t71, t59, t49 * t15 + t50 * t16, t69 * t8 + (-t16 * t41 + t44 * t46) * pkin(4) - t51 * t7, t15, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t45 * pkin(1) + t27 * pkin(2) + t8 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t59 - t15 * t39 + t16 * t33 + t74 * t46 + t69 * t7, -t65, t49 * t11 + t50 * t12, t69 * t4 + (-t12 * t41 + t21 * t44) * pkin(4) + t51 * t71, t11, -t76 * r_i_i_C(1) + t75 * r_i_i_C(2); 0, t64, t49 * t18 + t50 * t19, t69 * t10 + (-t19 * t41 + t24 * t44) * pkin(4) + t51 * (-t19 * t34 + t24 * t35), t18, (-t10 * t40 + t18 * t43) * r_i_i_C(1) + (-t10 * t43 - t18 * t40) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end