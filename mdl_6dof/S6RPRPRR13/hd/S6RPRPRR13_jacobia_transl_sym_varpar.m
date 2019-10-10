% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR13_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR13_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR13_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
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
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
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
	t12 = [(-t4 * t9 - t10) * r_i_i_C(1) + (-t4 * t8 + t11) * r_i_i_C(2) - t5 * pkin(1) + t6 * t7, t5 * t2, 0, 0, 0, 0; (-t4 * t11 + t8) * r_i_i_C(1) + (-t4 * t10 - t9) * r_i_i_C(2) + t6 * pkin(1) + t5 * t7, -t6 * t2, 0, 0, 0, 0; 0, t4, 0, 0, 0, 0;];
	Ja_transl = t12;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.28s
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
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (118->34), mult. (317->57), div. (0->0), fcn. (413->10), ass. (0->34)
	t44 = r_i_i_C(1) + pkin(9);
	t21 = cos(pkin(6));
	t16 = sin(pkin(12));
	t25 = cos(qJ(1));
	t34 = t25 * t16;
	t19 = cos(pkin(12));
	t23 = sin(qJ(1));
	t35 = t23 * t19;
	t11 = t21 * t34 + t35;
	t22 = sin(qJ(3));
	t24 = cos(qJ(3));
	t32 = t25 * t19;
	t37 = t23 * t16;
	t10 = -t21 * t32 + t37;
	t17 = sin(pkin(7));
	t20 = cos(pkin(7));
	t18 = sin(pkin(6));
	t33 = t25 * t18;
	t27 = t10 * t20 + t17 * t33;
	t43 = -t11 * t24 + t27 * t22;
	t1 = t11 * t22 + t27 * t24;
	t42 = -r_i_i_C(2) + pkin(3);
	t39 = t17 * t21;
	t38 = t20 * t24;
	t36 = t23 * t18;
	t31 = r_i_i_C(3) + qJ(4);
	t30 = t18 * qJ(2);
	t29 = t17 * t36;
	t13 = -t21 * t37 + t32;
	t12 = -t21 * t35 - t34;
	t7 = -t24 * t39 + (t16 * t22 - t19 * t38) * t18;
	t6 = t13 * t24 + (t12 * t20 + t29) * t22;
	t5 = -t12 * t38 + t13 * t22 - t24 * t29;
	t2 = [-t11 * pkin(2) - t23 * pkin(1) + t25 * t30 + t42 * t43 - t31 * t1 + t44 * (-t10 * t17 + t20 * t33), t36, t31 * t6 - t42 * t5, t5, 0, 0; t25 * pkin(1) + t13 * pkin(2) + t23 * t30 + t31 * t5 + t42 * t6 + t44 * (-t12 * t17 + t20 * t36), -t33, -t42 * t1 - t31 * t43, t1, 0, 0; 0, t21, t31 * (t22 * t39 + (t19 * t20 * t22 + t16 * t24) * t18) - t42 * t7, t7, 0, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (215->51), mult. (585->86), div. (0->0), fcn. (765->12), ass. (0->44)
	t55 = pkin(4) + pkin(9);
	t30 = cos(pkin(6));
	t25 = sin(pkin(12));
	t36 = cos(qJ(1));
	t45 = t36 * t25;
	t28 = cos(pkin(12));
	t33 = sin(qJ(1));
	t46 = t33 * t28;
	t18 = t30 * t45 + t46;
	t32 = sin(qJ(3));
	t35 = cos(qJ(3));
	t43 = t36 * t28;
	t48 = t33 * t25;
	t17 = -t30 * t43 + t48;
	t29 = cos(pkin(7));
	t26 = sin(pkin(7));
	t27 = sin(pkin(6));
	t44 = t36 * t27;
	t39 = t26 * t44;
	t37 = t17 * t29 + t39;
	t54 = -t18 * t35 + t32 * t37;
	t31 = sin(qJ(5));
	t34 = cos(qJ(5));
	t38 = t31 * r_i_i_C(1) + t34 * r_i_i_C(2) + qJ(4);
	t53 = t18 * t32;
	t51 = t26 * t30;
	t50 = t27 * t28;
	t49 = t29 * t35;
	t47 = t33 * t27;
	t42 = t27 * qJ(2);
	t41 = r_i_i_C(3) + pkin(10) + pkin(3);
	t40 = t26 * t47;
	t11 = -t17 * t26 + t29 * t44;
	t19 = -t30 * t46 - t45;
	t13 = -t19 * t26 + t29 * t47;
	t20 = -t30 * t48 + t43;
	t16 = -t26 * t50 + t30 * t29;
	t9 = t27 * t25 * t32 - t35 * t51 - t49 * t50;
	t8 = t20 * t35 + (t19 * t29 + t40) * t32;
	t7 = -t19 * t49 + t20 * t32 - t35 * t40;
	t3 = t17 * t49 + t35 * t39 + t53;
	t2 = t13 * t34 + t7 * t31;
	t1 = -t13 * t31 + t7 * t34;
	t4 = [-t18 * pkin(2) - t33 * pkin(1) + t36 * t42 + t41 * t54 + t38 * (-t35 * t37 - t53) + (t34 * r_i_i_C(1) - t31 * r_i_i_C(2) + t55) * t11, t47, t38 * t8 - t41 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t36 * pkin(1) + t20 * pkin(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(4) + t55 * t13 + t33 * t42 + t41 * t8, -t44, -t41 * t3 - t38 * t54, t3, (t11 * t31 + t3 * t34) * r_i_i_C(1) + (t11 * t34 - t3 * t31) * r_i_i_C(2), 0; 0, t30, -t41 * t9 + t38 * (t32 * t51 + (t28 * t29 * t32 + t25 * t35) * t27), t9, (-t16 * t31 + t9 * t34) * r_i_i_C(1) + (-t16 * t34 - t9 * t31) * r_i_i_C(2), 0;];
	Ja_transl = t4;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.39s
	% Computational Cost: add. (453->66), mult. (1247->111), div. (0->0), fcn. (1643->14), ass. (0->52)
	t67 = pkin(4) + pkin(9);
	t37 = cos(pkin(6));
	t35 = cos(pkin(12));
	t44 = cos(qJ(1));
	t56 = t44 * t35;
	t32 = sin(pkin(12));
	t41 = sin(qJ(1));
	t61 = t41 * t32;
	t26 = -t37 * t56 + t61;
	t58 = t44 * t32;
	t59 = t41 * t35;
	t27 = t37 * t58 + t59;
	t33 = sin(pkin(7));
	t36 = cos(pkin(7));
	t40 = sin(qJ(3));
	t34 = sin(pkin(6));
	t57 = t44 * t34;
	t63 = cos(qJ(3));
	t66 = (t26 * t36 + t33 * t57) * t40 - t27 * t63;
	t39 = sin(qJ(5));
	t43 = cos(qJ(5));
	t38 = sin(qJ(6));
	t42 = cos(qJ(6));
	t50 = t42 * r_i_i_C(1) - t38 * r_i_i_C(2) + pkin(5);
	t64 = r_i_i_C(3) + pkin(11);
	t45 = t50 * t39 - t64 * t43 + qJ(4);
	t65 = pkin(10) + pkin(3);
	t62 = t34 * t35;
	t60 = t41 * t34;
	t55 = t34 * qJ(2);
	t53 = t33 * t63;
	t52 = t36 * t63;
	t51 = t34 * t53;
	t20 = -t26 * t33 + t36 * t57;
	t48 = t37 * t59 + t58;
	t47 = t38 * r_i_i_C(1) + t42 * r_i_i_C(2) + t65;
	t46 = t48 * t36;
	t22 = t48 * t33 + t36 * t60;
	t12 = t26 * t52 + t27 * t40 + t44 * t51;
	t28 = -t37 * t61 + t56;
	t25 = -t33 * t62 + t37 * t36;
	t19 = t37 * t33 * t40 + (t35 * t36 * t40 + t63 * t32) * t34;
	t18 = t34 * t32 * t40 - t37 * t53 - t52 * t62;
	t17 = t28 * t63 + (t33 * t60 - t46) * t40;
	t16 = t28 * t40 - t41 * t51 + t63 * t46;
	t11 = t18 * t39 + t25 * t43;
	t8 = t16 * t39 + t22 * t43;
	t7 = -t16 * t43 + t22 * t39;
	t6 = t12 * t39 - t20 * t43;
	t2 = t17 * t38 + t8 * t42;
	t1 = t17 * t42 - t8 * t38;
	t3 = [t44 * t55 - t41 * pkin(1) - t27 * pkin(2) + t47 * t66 - t45 * t12 + (t64 * t39 + t50 * t43 + t67) * t20, t60, -t16 * t47 + t45 * t17, t16, -t50 * t7 + t64 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t44 * pkin(1) + t28 * pkin(2) + t8 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * qJ(4) + t65 * t17 + t67 * t22 + t41 * t55 + t64 * t7, -t57, -t12 * t47 - t45 * t66, t12, t64 * t6 + t50 * (t12 * t43 + t20 * t39), (-t6 * t38 - t42 * t66) * r_i_i_C(1) + (t38 * t66 - t6 * t42) * r_i_i_C(2); 0, t37, -t18 * t47 + t45 * t19, t18, t64 * t11 + t50 * (t18 * t43 - t25 * t39), (-t11 * t38 + t19 * t42) * r_i_i_C(1) + (-t11 * t42 - t19 * t38) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end