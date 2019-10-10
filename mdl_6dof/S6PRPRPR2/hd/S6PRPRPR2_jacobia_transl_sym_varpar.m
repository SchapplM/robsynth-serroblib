% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobia_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:50
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:50
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:50
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (7->7), mult. (20->17), div. (0->0), fcn. (24->6), ass. (0->8)
	t4 = cos(pkin(6));
	t5 = sin(qJ(2));
	t8 = t4 * t5;
	t6 = cos(qJ(2));
	t7 = t4 * t6;
	t3 = cos(pkin(10));
	t1 = sin(pkin(10));
	t2 = [0, (-t1 * t7 - t3 * t5) * r_i_i_C(1) + (t1 * t8 - t3 * t6) * r_i_i_C(2), 0, 0, 0, 0; 0, (-t1 * t5 + t3 * t7) * r_i_i_C(1) + (-t1 * t6 - t3 * t8) * r_i_i_C(2), 0, 0, 0, 0; 1, (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * sin(pkin(6)), 0, 0, 0, 0;];
	Ja_transl = t2;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:50
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (23->15), mult. (60->31), div. (0->0), fcn. (77->8), ass. (0->14)
	t10 = cos(pkin(6));
	t12 = cos(qJ(2));
	t13 = t10 * t12;
	t11 = sin(qJ(2));
	t5 = sin(pkin(11));
	t8 = cos(pkin(11));
	t4 = -t11 * t8 - t12 * t5;
	t3 = t11 * t5 - t12 * t8;
	t9 = cos(pkin(10));
	t7 = sin(pkin(6));
	t6 = sin(pkin(10));
	t2 = t4 * t10;
	t1 = t3 * t10;
	t14 = [0, (t6 * t1 + t9 * t4) * r_i_i_C(1) + (-t6 * t2 + t9 * t3) * r_i_i_C(2) + (-t9 * t11 - t6 * t13) * pkin(2), t6 * t7, 0, 0, 0; 0, (-t9 * t1 + t6 * t4) * r_i_i_C(1) + (t9 * t2 + t6 * t3) * r_i_i_C(2) + (-t6 * t11 + t9 * t13) * pkin(2), -t9 * t7, 0, 0, 0; 1, (t12 * pkin(2) - t3 * r_i_i_C(1) + t4 * r_i_i_C(2)) * t7, t10, 0, 0, 0;];
	Ja_transl = t14;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:50
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (76->27), mult. (197->52), div. (0->0), fcn. (255->10), ass. (0->23)
	t30 = pkin(8) + r_i_i_C(3);
	t15 = sin(pkin(10));
	t16 = sin(pkin(6));
	t29 = t15 * t16;
	t18 = cos(pkin(10));
	t28 = t18 * t16;
	t19 = cos(pkin(6));
	t23 = cos(qJ(2));
	t27 = t19 * t23;
	t14 = sin(pkin(11));
	t17 = cos(pkin(11));
	t21 = sin(qJ(2));
	t25 = t14 * t23 + t21 * t17;
	t10 = t25 * t19;
	t11 = t14 * t21 - t23 * t17;
	t3 = t10 * t18 - t11 * t15;
	t26 = t10 * t15 + t11 * t18;
	t20 = sin(qJ(4));
	t22 = cos(qJ(4));
	t24 = r_i_i_C(1) * t22 - r_i_i_C(2) * t20 + pkin(3);
	t9 = t11 * t19;
	t8 = t25 * t16;
	t1 = [0, -t30 * t26 + (-t15 * t27 - t18 * t21) * pkin(2) + t24 * (t15 * t9 - t18 * t25), t29, (t20 * t26 + t22 * t29) * r_i_i_C(1) + (-t20 * t29 + t22 * t26) * r_i_i_C(2), 0, 0; 0, t30 * t3 + (-t15 * t21 + t18 * t27) * pkin(2) + t24 * (-t15 * t25 - t18 * t9), -t28, (-t20 * t3 - t22 * t28) * r_i_i_C(1) + (t20 * t28 - t22 * t3) * r_i_i_C(2), 0, 0; 1, t30 * t8 + (pkin(2) * t23 - t11 * t24) * t16, t19, (t19 * t22 - t20 * t8) * r_i_i_C(1) + (-t19 * t20 - t22 * t8) * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (165->31), mult. (431->56), div. (0->0), fcn. (566->12), ass. (0->30)
	t25 = sin(pkin(10));
	t26 = sin(pkin(6));
	t43 = t25 * t26;
	t29 = cos(pkin(10));
	t42 = t29 * t26;
	t30 = cos(pkin(6));
	t34 = cos(qJ(2));
	t41 = t30 * t34;
	t40 = r_i_i_C(3) + qJ(5);
	t24 = sin(pkin(11));
	t28 = cos(pkin(11));
	t32 = sin(qJ(2));
	t38 = t34 * t24 + t32 * t28;
	t17 = t38 * t30;
	t18 = t32 * t24 - t34 * t28;
	t7 = t29 * t17 - t25 * t18;
	t39 = t25 * t17 + t29 * t18;
	t23 = sin(pkin(12));
	t27 = cos(pkin(12));
	t37 = r_i_i_C(1) * t27 - r_i_i_C(2) * t23 + pkin(4);
	t36 = t23 * r_i_i_C(1) + t27 * r_i_i_C(2) + pkin(8);
	t31 = sin(qJ(4));
	t33 = cos(qJ(4));
	t35 = t40 * t31 + t37 * t33 + pkin(3);
	t16 = t18 * t30;
	t15 = t38 * t26;
	t11 = t15 * t31 - t30 * t33;
	t3 = -t31 * t39 - t33 * t43;
	t1 = t7 * t31 + t33 * t42;
	t2 = [0, (-t25 * t41 - t29 * t32) * pkin(2) - t36 * t39 + t35 * (t25 * t16 - t29 * t38), t43, t40 * (t31 * t43 - t33 * t39) - t37 * t3, t3, 0; 0, (-t25 * t32 + t29 * t41) * pkin(2) + t36 * t7 + t35 * (-t29 * t16 - t25 * t38), -t42, t40 * (-t31 * t42 + t7 * t33) - t37 * t1, t1, 0; 1, t36 * t15 + (pkin(2) * t34 - t18 * t35) * t26, t30, t40 * (t15 * t33 + t30 * t31) - t37 * t11, t11, 0;];
	Ja_transl = t2;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:50
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (260->43), mult. (577->77), div. (0->0), fcn. (759->14), ass. (0->37)
	t29 = sin(pkin(11));
	t36 = sin(qJ(2));
	t38 = cos(qJ(2));
	t46 = cos(pkin(11));
	t42 = -t36 * t29 + t38 * t46;
	t35 = sin(qJ(4));
	t37 = cos(qJ(4));
	t27 = pkin(12) + qJ(6);
	t25 = sin(t27);
	t26 = cos(t27);
	t43 = t26 * r_i_i_C(1) - t25 * r_i_i_C(2) + cos(pkin(12)) * pkin(5) + pkin(4);
	t51 = r_i_i_C(3) + pkin(9) + qJ(5);
	t39 = t51 * t35 + t43 * t37 + pkin(3);
	t30 = sin(pkin(10));
	t31 = sin(pkin(6));
	t50 = t30 * t31;
	t32 = cos(pkin(10));
	t49 = t32 * t31;
	t33 = cos(pkin(6));
	t48 = t33 * t38;
	t19 = -t38 * t29 - t36 * t46;
	t17 = t19 * t33;
	t7 = -t32 * t17 + t30 * t42;
	t44 = -t30 * t17 - t32 * t42;
	t41 = sin(pkin(12)) * pkin(5) + t25 * r_i_i_C(1) + t26 * r_i_i_C(2) + pkin(8);
	t40 = t42 * t33;
	t16 = t19 * t31;
	t15 = t42 * t31;
	t12 = -t16 * t37 + t33 * t35;
	t11 = -t16 * t35 - t33 * t37;
	t9 = t32 * t19 - t30 * t40;
	t6 = t30 * t19 + t32 * t40;
	t4 = t35 * t50 - t37 * t44;
	t3 = -t35 * t44 - t37 * t50;
	t2 = -t35 * t49 + t7 * t37;
	t1 = t7 * t35 + t37 * t49;
	t5 = [0, (-t30 * t48 - t32 * t36) * pkin(2) - t41 * t44 + t39 * t9, t50, -t43 * t3 + t51 * t4, t3, (-t4 * t25 - t9 * t26) * r_i_i_C(1) + (t9 * t25 - t4 * t26) * r_i_i_C(2); 0, (-t30 * t36 + t32 * t48) * pkin(2) + t41 * t7 + t39 * t6, -t49, -t43 * t1 + t51 * t2, t1, (-t2 * t25 - t6 * t26) * r_i_i_C(1) + (-t2 * t26 + t6 * t25) * r_i_i_C(2); 1, t31 * t38 * pkin(2) + t39 * t15 - t41 * t16, t33, -t43 * t11 + t51 * t12, t11, (-t12 * t25 - t15 * t26) * r_i_i_C(1) + (-t12 * t26 + t15 * t25) * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end