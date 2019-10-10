% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:40
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:40
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(11));
	t1 = sin(pkin(11));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:40
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (35->19), mult. (92->37), div. (0->0), fcn. (112->8), ass. (0->19)
	t6 = sin(pkin(6));
	t9 = sin(qJ(3));
	t20 = t6 * t9;
	t19 = pkin(8) + r_i_i_C(3);
	t11 = cos(qJ(3));
	t18 = t11 * t6;
	t10 = sin(qJ(2));
	t5 = sin(pkin(11));
	t17 = t5 * t10;
	t12 = cos(qJ(2));
	t16 = t5 * t12;
	t7 = cos(pkin(11));
	t15 = t7 * t10;
	t14 = t7 * t12;
	t13 = t11 * r_i_i_C(1) - t9 * r_i_i_C(2) + pkin(2);
	t8 = cos(pkin(6));
	t4 = -t8 * t17 + t14;
	t2 = t8 * t15 + t16;
	t1 = [0, t19 * t4 + t13 * (-t8 * t16 - t15), (t5 * t18 - t4 * t9) * r_i_i_C(1) + (-t11 * t4 - t5 * t20) * r_i_i_C(2), 0, 0, 0; 0, t19 * t2 + t13 * (t8 * t14 - t17), (-t7 * t18 - t2 * t9) * r_i_i_C(1) + (-t11 * t2 + t7 * t20) * r_i_i_C(2), 0, 0, 0; 1, (t19 * t10 + t13 * t12) * t6, (-t10 * t20 + t11 * t8) * r_i_i_C(1) + (-t10 * t18 - t8 * t9) * r_i_i_C(2), 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:40
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (102->32), mult. (269->63), div. (0->0), fcn. (338->10), ass. (0->28)
	t15 = sin(qJ(3));
	t18 = cos(qJ(3));
	t14 = sin(qJ(4));
	t17 = cos(qJ(4));
	t22 = t17 * r_i_i_C(1) - t14 * r_i_i_C(2) + pkin(3);
	t31 = pkin(9) + r_i_i_C(3);
	t32 = t31 * t15 + t22 * t18 + pkin(2);
	t13 = sin(pkin(6));
	t30 = t13 * t15;
	t29 = t13 * t18;
	t19 = cos(qJ(2));
	t28 = t13 * t19;
	t27 = cos(pkin(6));
	t26 = cos(pkin(11));
	t12 = sin(pkin(11));
	t25 = t12 * t27;
	t24 = t13 * t26;
	t23 = t27 * t26;
	t21 = t14 * r_i_i_C(1) + t17 * r_i_i_C(2) + pkin(8);
	t16 = sin(qJ(2));
	t10 = t27 * t15 + t16 * t29;
	t8 = -t16 * t25 + t26 * t19;
	t7 = t26 * t16 + t19 * t25;
	t6 = t12 * t19 + t16 * t23;
	t5 = t12 * t16 - t19 * t23;
	t4 = t12 * t30 + t8 * t18;
	t2 = -t15 * t24 + t6 * t18;
	t1 = [0, t21 * t8 - t32 * t7, t31 * t4 + t22 * (t12 * t29 - t8 * t15), (-t4 * t14 + t7 * t17) * r_i_i_C(1) + (-t7 * t14 - t4 * t17) * r_i_i_C(2), 0, 0; 0, t21 * t6 - t32 * t5, t31 * t2 + t22 * (-t6 * t15 - t18 * t24), (-t2 * t14 + t5 * t17) * r_i_i_C(1) + (-t5 * t14 - t2 * t17) * r_i_i_C(2), 0, 0; 1, (t21 * t16 + t32 * t19) * t13, t31 * t10 + t22 * (-t16 * t30 + t27 * t18), (-t10 * t14 - t17 * t28) * r_i_i_C(1) + (-t10 * t17 + t14 * t28) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:40
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (197->41), mult. (384->74), div. (0->0), fcn. (482->12), ass. (0->34)
	t25 = sin(qJ(3));
	t28 = cos(qJ(3));
	t21 = qJ(4) + qJ(5);
	t19 = sin(t21);
	t20 = cos(t21);
	t27 = cos(qJ(4));
	t33 = t27 * pkin(4) + r_i_i_C(1) * t20 - r_i_i_C(2) * t19 + pkin(3);
	t42 = r_i_i_C(3) + pkin(10) + pkin(9);
	t46 = t42 * t25 + t33 * t28 + pkin(2);
	t22 = sin(pkin(11));
	t26 = sin(qJ(2));
	t29 = cos(qJ(2));
	t37 = cos(pkin(11));
	t38 = cos(pkin(6));
	t34 = t38 * t37;
	t11 = t22 * t26 - t29 * t34;
	t12 = t22 * t29 + t26 * t34;
	t23 = sin(pkin(6));
	t35 = t23 * t37;
	t8 = t12 * t28 - t25 * t35;
	t45 = (t11 * t20 - t8 * t19) * r_i_i_C(1) + (-t11 * t19 - t8 * t20) * r_i_i_C(2);
	t36 = t22 * t38;
	t14 = -t26 * t36 + t37 * t29;
	t41 = t23 * t25;
	t10 = t14 * t28 + t22 * t41;
	t13 = t37 * t26 + t29 * t36;
	t44 = (-t10 * t19 + t13 * t20) * r_i_i_C(1) + (-t10 * t20 - t13 * t19) * r_i_i_C(2);
	t40 = t23 * t28;
	t16 = t38 * t25 + t26 * t40;
	t39 = t23 * t29;
	t43 = (-t16 * t19 - t20 * t39) * r_i_i_C(1) + (-t16 * t20 + t19 * t39) * r_i_i_C(2);
	t24 = sin(qJ(4));
	t32 = t24 * pkin(4) + t19 * r_i_i_C(1) + t20 * r_i_i_C(2) + pkin(8);
	t1 = [0, -t13 * t46 + t32 * t14, t42 * t10 + t33 * (-t14 * t25 + t22 * t40), (-t10 * t24 + t13 * t27) * pkin(4) + t44, t44, 0; 0, -t11 * t46 + t32 * t12, t42 * t8 + t33 * (-t12 * t25 - t28 * t35), (t11 * t27 - t24 * t8) * pkin(4) + t45, t45, 0; 1, (t32 * t26 + t46 * t29) * t23, t42 * t16 + t33 * (-t26 * t41 + t38 * t28), (-t16 * t24 - t27 * t39) * pkin(4) + t43, t43, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:07:40
	% EndTime: 2019-10-09 23:07:40
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (343->53), mult. (630->93), div. (0->0), fcn. (805->12), ass. (0->43)
	t69 = pkin(5) + r_i_i_C(1);
	t61 = r_i_i_C(3) + qJ(6);
	t45 = cos(qJ(4));
	t36 = t45 * pkin(4) + pkin(3);
	t43 = sin(qJ(3));
	t46 = cos(qJ(3));
	t68 = r_i_i_C(2) + pkin(10) + pkin(9);
	t70 = t36 * t46 + t68 * t43 + pkin(2);
	t39 = qJ(4) + qJ(5);
	t37 = sin(t39);
	t67 = t37 * t46;
	t38 = cos(t39);
	t66 = t38 * t46;
	t41 = sin(pkin(6));
	t65 = t41 * t43;
	t64 = t41 * t46;
	t47 = cos(qJ(2));
	t63 = t41 * t47;
	t62 = t46 * t47;
	t60 = cos(pkin(6));
	t59 = cos(pkin(11));
	t42 = sin(qJ(4));
	t58 = pkin(4) * t42 + pkin(8);
	t40 = sin(pkin(11));
	t56 = t40 * t60;
	t55 = t41 * t59;
	t44 = sin(qJ(2));
	t52 = t60 * t59;
	t29 = t40 * t47 + t44 * t52;
	t21 = t29 * t46 - t43 * t55;
	t28 = t40 * t44 - t47 * t52;
	t7 = t21 * t37 - t28 * t38;
	t54 = t61 * (t21 * t38 + t28 * t37) - t69 * t7;
	t31 = -t44 * t56 + t59 * t47;
	t23 = t31 * t46 + t40 * t65;
	t30 = t59 * t44 + t47 * t56;
	t9 = t23 * t37 - t30 * t38;
	t53 = -t69 * t9 + t61 * (t23 * t38 + t30 * t37);
	t33 = t60 * t43 + t44 * t64;
	t18 = t33 * t37 + t38 * t63;
	t51 = t61 * (t33 * t38 - t37 * t63) - t69 * t18;
	t49 = t61 * t37 + t69 * t38 + t36;
	t1 = [0, t58 * t31 + t69 * (-t30 * t66 + t31 * t37) + t61 * (-t30 * t67 - t31 * t38) - t70 * t30, t68 * t23 + t49 * (-t31 * t43 + t40 * t64), (-t23 * t42 + t30 * t45) * pkin(4) + t53, t53, t9; 0, t58 * t29 + t69 * (-t28 * t66 + t29 * t37) + t61 * (-t28 * t67 - t29 * t38) - t70 * t28, t68 * t21 + t49 * (-t29 * t43 - t46 * t55), (-t21 * t42 + t28 * t45) * pkin(4) + t54, t54, t7; 1, (t69 * (t37 * t44 + t38 * t62) + t61 * (t37 * t62 - t38 * t44) + t58 * t44 + t70 * t47) * t41, t68 * t33 + t49 * (-t44 * t65 + t60 * t46), (-t33 * t42 - t45 * t63) * pkin(4) + t51, t51, t18;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end