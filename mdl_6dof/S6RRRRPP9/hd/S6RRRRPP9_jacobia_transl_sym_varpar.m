% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:33
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP9_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP9_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP9_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:10
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:10
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
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:10
	% DurationCPUTime: 0.14s
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
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:10
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (150->46), mult. (381->80), div. (0->0), fcn. (482->10), ass. (0->34)
	t21 = sin(qJ(2));
	t22 = sin(qJ(1));
	t25 = cos(qJ(2));
	t26 = cos(qJ(1));
	t33 = cos(pkin(6));
	t31 = t26 * t33;
	t11 = t22 * t21 - t25 * t31;
	t19 = sin(qJ(4));
	t23 = cos(qJ(4));
	t12 = t21 * t31 + t22 * t25;
	t20 = sin(qJ(3));
	t24 = cos(qJ(3));
	t18 = sin(pkin(6));
	t34 = t18 * t26;
	t4 = t12 * t24 - t20 * t34;
	t43 = -t11 * t23 + t4 * t19;
	t42 = -t11 * t19 - t4 * t23;
	t30 = t23 * r_i_i_C(1) - t19 * r_i_i_C(2) + pkin(3);
	t40 = pkin(10) + r_i_i_C(3);
	t41 = t40 * t20 + t30 * t24 + pkin(2);
	t37 = t18 * t22;
	t36 = t18 * t24;
	t35 = t18 * t25;
	t32 = t22 * t33;
	t29 = t19 * r_i_i_C(1) + t23 * r_i_i_C(2) + pkin(9);
	t28 = -t12 * t20 - t24 * t34;
	t14 = -t21 * t32 + t26 * t25;
	t13 = t26 * t21 + t25 * t32;
	t10 = t33 * t20 + t21 * t36;
	t8 = t14 * t24 + t20 * t37;
	t7 = t14 * t20 - t22 * t36;
	t2 = t13 * t19 + t8 * t23;
	t1 = t13 * t23 - t8 * t19;
	t3 = [-t22 * pkin(1) - t12 * pkin(2) - t4 * pkin(3) + pkin(8) * t34 - t11 * pkin(9) + t42 * r_i_i_C(1) + t43 * r_i_i_C(2) + t40 * t28, -t13 * t41 + t29 * t14, -t30 * t7 + t40 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t26 * pkin(1) + t14 * pkin(2) + t8 * pkin(3) + pkin(8) * t37 + t13 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t40 * t7, -t11 * t41 + t29 * t12, t30 * t28 + t40 * t4, -t43 * r_i_i_C(1) + t42 * r_i_i_C(2), 0, 0; 0, (t29 * t21 + t41 * t25) * t18, t40 * t10 + t30 * (-t18 * t21 * t20 + t33 * t24), (-t10 * t19 - t23 * t35) * r_i_i_C(1) + (-t10 * t23 + t19 * t35) * r_i_i_C(2), 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:10
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (242->58), mult. (613->99), div. (0->0), fcn. (787->10), ass. (0->39)
	t34 = sin(qJ(2));
	t35 = sin(qJ(1));
	t38 = cos(qJ(2));
	t39 = cos(qJ(1));
	t46 = cos(pkin(6));
	t43 = t39 * t46;
	t25 = t34 * t43 + t35 * t38;
	t33 = sin(qJ(3));
	t37 = cos(qJ(3));
	t31 = sin(pkin(6));
	t51 = t31 * t39;
	t14 = t25 * t37 - t33 * t51;
	t24 = t35 * t34 - t38 * t43;
	t32 = sin(qJ(4));
	t36 = cos(qJ(4));
	t1 = t14 * t32 - t24 * t36;
	t60 = t14 * t36 + t24 * t32;
	t57 = pkin(10) + r_i_i_C(1);
	t59 = pkin(3) * t37 + t57 * t33 + pkin(2);
	t47 = r_i_i_C(3) + qJ(5);
	t58 = pkin(4) - r_i_i_C(2);
	t40 = t47 * t32 + t58 * t36 + pkin(3);
	t54 = t31 * t35;
	t53 = t31 * t37;
	t52 = t31 * t38;
	t50 = t32 * t37;
	t49 = t36 * t37;
	t48 = t37 * t38;
	t44 = t35 * t46;
	t42 = -t25 * t33 - t37 * t51;
	t27 = -t34 * t44 + t39 * t38;
	t26 = t39 * t34 + t38 * t44;
	t23 = t46 * t33 + t34 * t53;
	t18 = t27 * t37 + t33 * t54;
	t17 = t27 * t33 - t35 * t53;
	t11 = t23 * t32 + t36 * t52;
	t6 = t18 * t36 + t26 * t32;
	t5 = t18 * t32 - t26 * t36;
	t2 = [-t35 * pkin(1) - t25 * pkin(2) - t14 * pkin(3) + pkin(8) * t51 - t24 * pkin(9) - t47 * t1 + t57 * t42 - t58 * t60, t27 * pkin(9) + t47 * (-t26 * t50 - t27 * t36) + t58 * (-t26 * t49 + t27 * t32) - t59 * t26, -t40 * t17 + t57 * t18, t47 * t6 - t58 * t5, t5, 0; t39 * pkin(1) + t27 * pkin(2) + t18 * pkin(3) + pkin(8) * t54 + t26 * pkin(9) + t57 * t17 + t47 * t5 + t58 * t6, t25 * pkin(9) + t58 * (-t24 * t49 + t25 * t32) + t47 * (-t24 * t50 - t25 * t36) - t59 * t24, t57 * t14 + t40 * t42, -t58 * t1 + t47 * t60, t1, 0; 0, (t58 * (t32 * t34 + t36 * t48) + t47 * (t32 * t48 - t34 * t36) + pkin(9) * t34 + t59 * t38) * t31, t57 * t23 + t40 * (-t31 * t34 * t33 + t46 * t37), t47 * (t23 * t36 - t32 * t52) - t58 * t11, t11, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:33:10
	% EndTime: 2019-10-10 12:33:11
	% DurationCPUTime: 0.33s
	% Computational Cost: add. (312->58), mult. (788->99), div. (0->0), fcn. (1016->10), ass. (0->40)
	t35 = sin(qJ(2));
	t36 = sin(qJ(1));
	t39 = cos(qJ(2));
	t49 = cos(pkin(6));
	t59 = cos(qJ(1));
	t43 = t49 * t59;
	t26 = t35 * t43 + t36 * t39;
	t34 = sin(qJ(3));
	t38 = cos(qJ(3));
	t32 = sin(pkin(6));
	t45 = t32 * t59;
	t14 = t26 * t38 - t34 * t45;
	t25 = t36 * t35 - t39 * t43;
	t33 = sin(qJ(4));
	t37 = cos(qJ(4));
	t1 = t14 * t33 - t25 * t37;
	t2 = t14 * t37 + t25 * t33;
	t48 = pkin(5) + pkin(10) + r_i_i_C(1);
	t60 = pkin(3) * t38 + t48 * t34 + pkin(2);
	t47 = pkin(4) + r_i_i_C(3) + qJ(6);
	t50 = r_i_i_C(2) + qJ(5);
	t40 = t50 * t33 + t47 * t37 + pkin(3);
	t56 = t32 * t36;
	t55 = t32 * t38;
	t54 = t32 * t39;
	t53 = t33 * t38;
	t52 = t37 * t38;
	t51 = t38 * t39;
	t44 = t36 * t49;
	t42 = -t26 * t34 - t38 * t45;
	t28 = -t35 * t44 + t59 * t39;
	t27 = t59 * t35 + t39 * t44;
	t24 = t49 * t34 + t35 * t55;
	t18 = t28 * t38 + t34 * t56;
	t17 = t28 * t34 - t36 * t55;
	t12 = t24 * t37 - t33 * t54;
	t11 = t24 * t33 + t37 * t54;
	t6 = t18 * t37 + t27 * t33;
	t5 = t18 * t33 - t27 * t37;
	t3 = [-t36 * pkin(1) - t26 * pkin(2) - t14 * pkin(3) + pkin(8) * t45 - t25 * pkin(9) - t50 * t1 - t47 * t2 + t48 * t42, t28 * pkin(9) + t50 * (-t27 * t53 - t28 * t37) + t47 * (-t27 * t52 + t28 * t33) - t60 * t27, -t40 * t17 + t48 * t18, -t47 * t5 + t50 * t6, t5, t6; t59 * pkin(1) + t28 * pkin(2) + t18 * pkin(3) + pkin(8) * t56 + t27 * pkin(9) + t48 * t17 + t47 * t6 + t50 * t5, t26 * pkin(9) + t50 * (-t25 * t53 - t26 * t37) + t47 * (-t25 * t52 + t26 * t33) - t60 * t25, t48 * t14 + t40 * t42, -t47 * t1 + t50 * t2, t1, t2; 0, (t50 * (t33 * t51 - t35 * t37) + t47 * (t33 * t35 + t37 * t51) + pkin(9) * t35 + t60 * t39) * t32, t48 * t24 + t40 * (-t32 * t35 * t34 + t49 * t38), -t47 * t11 + t50 * t12, t11, t12;];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,6);
end