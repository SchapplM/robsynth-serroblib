% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR14
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRRR14_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR14_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR14_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:36
	% EndTime: 2019-12-29 18:05:37
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:37
	% EndTime: 2019-12-29 18:05:37
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
	% StartTime: 2019-12-29 18:05:36
	% EndTime: 2019-12-29 18:05:37
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (13->12), mult. (28->19), div. (0->0), fcn. (35->6), ass. (0->12)
	t1 = sin(pkin(11));
	t5 = sin(qJ(1));
	t11 = t5 * t1;
	t3 = cos(pkin(11));
	t10 = t5 * t3;
	t6 = cos(qJ(1));
	t9 = t6 * t1;
	t8 = t6 * t3;
	t2 = sin(pkin(5));
	t7 = t2 * (r_i_i_C(3) + qJ(2));
	t4 = cos(pkin(5));
	t12 = [(-t4 * t9 - t10) * r_i_i_C(1) + (-t4 * t8 + t11) * r_i_i_C(2) - t5 * pkin(1) + t6 * t7, t5 * t2, 0, 0, 0; (-t4 * t11 + t8) * r_i_i_C(1) + (-t4 * t10 - t9) * r_i_i_C(2) + t6 * pkin(1) + t5 * t7, -t6 * t2, 0, 0, 0; 0, t4, 0, 0, 0;];
	Ja_transl = t12;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:37
	% EndTime: 2019-12-29 18:05:37
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (62->32), mult. (166->58), div. (0->0), fcn. (213->10), ass. (0->29)
	t34 = r_i_i_C(3) + pkin(8);
	t13 = cos(pkin(6));
	t15 = sin(qJ(3));
	t17 = cos(qJ(3));
	t10 = sin(pkin(6));
	t11 = sin(pkin(5));
	t18 = cos(qJ(1));
	t26 = t18 * t11;
	t24 = t10 * t26;
	t12 = cos(pkin(11));
	t14 = cos(pkin(5));
	t29 = t14 * t18;
	t16 = sin(qJ(1));
	t9 = sin(pkin(11));
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
	t7 = [t19 * r_i_i_C(1) + t33 * r_i_i_C(2) - t4 * pkin(2) - t16 * pkin(1) + t18 * t25 + t34 * (-t3 * t10 + t13 * t26), t28, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0; t18 * pkin(1) + t6 * pkin(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t25 + t34 * (-t10 * t5 + t13 * t28), -t26, -t33 * r_i_i_C(1) + t19 * r_i_i_C(2), 0, 0; 0, t14, (r_i_i_C(1) * t17 - r_i_i_C(2) * t15) * t14 * t10 + ((t12 * t13 * t17 - t15 * t9) * r_i_i_C(1) + (-t12 * t30 - t17 * t9) * r_i_i_C(2)) * t11, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:37
	% EndTime: 2019-12-29 18:05:37
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (176->48), mult. (480->84), div. (0->0), fcn. (625->12), ass. (0->42)
	t27 = cos(pkin(5));
	t25 = cos(pkin(11));
	t33 = cos(qJ(1));
	t41 = t33 * t25;
	t22 = sin(pkin(11));
	t30 = sin(qJ(1));
	t46 = t30 * t22;
	t16 = -t27 * t41 + t46;
	t23 = sin(pkin(6));
	t26 = cos(pkin(6));
	t24 = sin(pkin(5));
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
	t49 = r_i_i_C(3) + pkin(9);
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
	t3 = [-t30 * pkin(1) - t17 * pkin(2) + t6 * pkin(3) + t11 * pkin(8) + t51 * r_i_i_C(1) + t52 * r_i_i_C(2) + t33 * t40 + t49 * t34, t45, -t37 * t7 + t49 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t33 * pkin(1) + t18 * pkin(2) + t8 * pkin(3) + t13 * pkin(8) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t30 * t40 + t49 * t7, -t42, t37 * t34 - t49 * t6, -t52 * r_i_i_C(1) + t51 * r_i_i_C(2), 0; 0, t27, t49 * t10 + t37 * (t32 * t48 + (t25 * t26 * t32 - t22 * t29) * t24), (-t10 * t28 + t15 * t31) * r_i_i_C(1) + (-t10 * t31 - t15 * t28) * r_i_i_C(2), 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:05:36
	% EndTime: 2019-12-29 18:05:37
	% DurationCPUTime: 0.67s
	% Computational Cost: add. (414->64), mult. (1142->110), div. (0->0), fcn. (1503->14), ass. (0->51)
	t39 = cos(qJ(1));
	t56 = sin(pkin(11));
	t60 = cos(pkin(5));
	t48 = t60 * t56;
	t58 = cos(pkin(11));
	t62 = sin(qJ(1));
	t26 = t39 * t48 + t62 * t58;
	t36 = sin(qJ(3));
	t63 = cos(qJ(3));
	t50 = t60 * t58;
	t25 = -t39 * t50 + t62 * t56;
	t33 = sin(pkin(5));
	t57 = sin(pkin(6));
	t52 = t33 * t57;
	t59 = cos(pkin(6));
	t67 = t25 * t59 + t39 * t52;
	t11 = t26 * t36 + t67 * t63;
	t34 = sin(qJ(5));
	t37 = cos(qJ(5));
	t12 = t26 * t63 - t67 * t36;
	t53 = t33 * t59;
	t21 = t25 * t57 - t39 * t53;
	t35 = sin(qJ(4));
	t38 = cos(qJ(4));
	t4 = t12 * t38 + t21 * t35;
	t71 = -t11 * t37 + t4 * t34;
	t70 = -t11 * t34 - t4 * t37;
	t66 = -t12 * t35 + t21 * t38;
	t42 = t39 * t56 + t62 * t50;
	t65 = t42 * t59 - t62 * t52;
	t64 = r_i_i_C(3) + pkin(10);
	t61 = t39 * t33;
	t55 = t62 * t33;
	t49 = t60 * t57;
	t47 = t59 * t58;
	t46 = t37 * r_i_i_C(1) - t34 * r_i_i_C(2) + pkin(4);
	t45 = t34 * r_i_i_C(1) + t37 * r_i_i_C(2) + pkin(9);
	t43 = -t64 * t35 - t46 * t38 - pkin(3);
	t40 = t42 * t57 + t62 * t53;
	t27 = t39 * t58 - t62 * t48;
	t24 = -t58 * t52 + t60 * t59;
	t19 = t36 * t49 + (t36 * t47 + t63 * t56) * t33;
	t18 = -t63 * t49 + (t36 * t56 - t47 * t63) * t33;
	t16 = t27 * t63 - t65 * t36;
	t15 = t27 * t36 + t65 * t63;
	t10 = t19 * t38 + t24 * t35;
	t8 = t16 * t38 + t40 * t35;
	t7 = t16 * t35 - t40 * t38;
	t2 = t15 * t34 + t8 * t37;
	t1 = t15 * t37 - t8 * t34;
	t3 = [-t62 * pkin(1) - t26 * pkin(2) - t12 * pkin(3) - t4 * pkin(4) - t21 * pkin(8) - t11 * pkin(9) + t70 * r_i_i_C(1) + t71 * r_i_i_C(2) + qJ(2) * t61 + t64 * t66, t55, t43 * t15 + t45 * t16, -t46 * t7 + t64 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t39 * pkin(1) + t27 * pkin(2) + t16 * pkin(3) + t8 * pkin(4) + t40 * pkin(8) + t15 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t55 + t64 * t7, -t61, t43 * t11 + t45 * t12, t64 * t4 + t46 * t66, -t71 * r_i_i_C(1) + t70 * r_i_i_C(2); 0, t60, t43 * t18 + t45 * t19, t64 * t10 + t46 * (-t19 * t35 + t24 * t38), (-t10 * t34 + t18 * t37) * r_i_i_C(1) + (-t10 * t37 - t18 * t34) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,5);
end