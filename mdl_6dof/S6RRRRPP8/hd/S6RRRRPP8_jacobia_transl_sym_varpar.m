% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRRPP8
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
% Datum: 2019-10-10 12:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPP8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:15
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:15
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
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:15
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
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.27s
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
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.41s
	% Computational Cost: add. (242->58), mult. (613->99), div. (0->0), fcn. (787->10), ass. (0->39)
	t35 = sin(qJ(2));
	t36 = sin(qJ(1));
	t39 = cos(qJ(2));
	t40 = cos(qJ(1));
	t47 = cos(pkin(6));
	t44 = t40 * t47;
	t26 = t35 * t44 + t36 * t39;
	t34 = sin(qJ(3));
	t38 = cos(qJ(3));
	t32 = sin(pkin(6));
	t52 = t32 * t40;
	t14 = t26 * t38 - t34 * t52;
	t25 = t36 * t35 - t39 * t44;
	t33 = sin(qJ(4));
	t37 = cos(qJ(4));
	t1 = t14 * t33 - t25 * t37;
	t61 = t14 * t37 + t25 * t33;
	t58 = pkin(10) + r_i_i_C(2);
	t60 = pkin(3) * t38 + t58 * t34 + pkin(2);
	t48 = r_i_i_C(3) + qJ(5);
	t59 = pkin(4) + r_i_i_C(1);
	t41 = t48 * t33 + t59 * t37 + pkin(3);
	t55 = t32 * t36;
	t54 = t32 * t38;
	t53 = t32 * t39;
	t51 = t33 * t38;
	t50 = t37 * t38;
	t49 = t38 * t39;
	t45 = t36 * t47;
	t43 = -t26 * t34 - t38 * t52;
	t28 = -t35 * t45 + t40 * t39;
	t27 = t40 * t35 + t39 * t45;
	t24 = t47 * t34 + t35 * t54;
	t18 = t28 * t38 + t34 * t55;
	t17 = t28 * t34 - t36 * t54;
	t11 = t24 * t33 + t37 * t53;
	t6 = t18 * t37 + t27 * t33;
	t5 = t18 * t33 - t27 * t37;
	t2 = [-t36 * pkin(1) - t26 * pkin(2) - t14 * pkin(3) + pkin(8) * t52 - t25 * pkin(9) - t48 * t1 + t58 * t43 - t59 * t61, t28 * pkin(9) + t48 * (-t27 * t51 - t28 * t37) + t59 * (-t27 * t50 + t28 * t33) - t60 * t27, -t41 * t17 + t58 * t18, t48 * t6 - t59 * t5, t5, 0; t40 * pkin(1) + t28 * pkin(2) + t18 * pkin(3) + pkin(8) * t55 + t27 * pkin(9) + t58 * t17 + t48 * t5 + t59 * t6, t26 * pkin(9) + t59 * (-t25 * t50 + t26 * t33) + t48 * (-t25 * t51 - t26 * t37) - t60 * t25, t58 * t14 + t41 * t43, -t59 * t1 + t48 * t61, t1, 0; 0, (t59 * (t33 * t35 + t37 * t49) + t48 * (t33 * t49 - t35 * t37) + pkin(9) * t35 + t60 * t39) * t32, t58 * t24 + t41 * (-t32 * t35 * t34 + t47 * t38), t48 * (t24 * t37 - t33 * t53) - t59 * t11, t11, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:31:15
	% EndTime: 2019-10-10 12:31:16
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (307->59), mult. (775->98), div. (0->0), fcn. (998->10), ass. (0->40)
	t36 = sin(qJ(2));
	t39 = cos(qJ(2));
	t49 = cos(pkin(6));
	t59 = cos(qJ(1));
	t43 = t49 * t59;
	t58 = sin(qJ(1));
	t26 = t36 * t43 + t58 * t39;
	t35 = sin(qJ(3));
	t38 = cos(qJ(3));
	t33 = sin(pkin(6));
	t45 = t33 * t59;
	t14 = t26 * t38 - t35 * t45;
	t25 = t58 * t36 - t39 * t43;
	t34 = sin(qJ(4));
	t37 = cos(qJ(4));
	t1 = t14 * t34 - t25 * t37;
	t61 = t14 * t37 + t25 * t34;
	t47 = pkin(10) - r_i_i_C(3) - qJ(6);
	t60 = pkin(3) * t38 + t47 * t35 + pkin(2);
	t48 = pkin(4) + pkin(5) + r_i_i_C(1);
	t50 = r_i_i_C(2) + qJ(5);
	t40 = t50 * t34 + t48 * t37 + pkin(3);
	t55 = t33 * t36;
	t54 = t33 * t39;
	t53 = t34 * t38;
	t52 = t37 * t38;
	t51 = t38 * t39;
	t44 = t33 * t58;
	t13 = -t26 * t35 - t38 * t45;
	t42 = t49 * t58;
	t28 = -t36 * t42 + t59 * t39;
	t27 = t59 * t36 + t39 * t42;
	t24 = t49 * t35 + t38 * t55;
	t23 = -t35 * t55 + t49 * t38;
	t18 = t28 * t38 + t35 * t44;
	t17 = t28 * t35 - t38 * t44;
	t11 = t24 * t34 + t37 * t54;
	t6 = t18 * t37 + t27 * t34;
	t5 = t18 * t34 - t27 * t37;
	t2 = [-t58 * pkin(1) - t26 * pkin(2) - t14 * pkin(3) + pkin(8) * t45 - t25 * pkin(9) - t50 * t1 + t47 * t13 - t48 * t61, t28 * pkin(9) + t50 * (-t27 * t53 - t28 * t37) + t48 * (-t27 * t52 + t28 * t34) - t60 * t27, -t40 * t17 + t47 * t18, -t48 * t5 + t50 * t6, t5, -t17; t59 * pkin(1) + t28 * pkin(2) + t18 * pkin(3) + pkin(8) * t44 + t27 * pkin(9) + t47 * t17 + t48 * t6 + t50 * t5, t26 * pkin(9) + t50 * (-t25 * t53 - t26 * t37) + t48 * (-t25 * t52 + t26 * t34) - t60 * t25, t40 * t13 + t47 * t14, -t48 * t1 + t50 * t61, t1, t13; 0, (t50 * (t34 * t51 - t36 * t37) + t48 * (t34 * t36 + t37 * t51) + pkin(9) * t36 + t60 * t39) * t33, t40 * t23 + t47 * t24, t50 * (t24 * t37 - t34 * t54) - t48 * t11, t11, t23;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,6);
end