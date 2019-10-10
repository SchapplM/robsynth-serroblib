% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:58
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:49
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:49
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
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:49
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t2 * t5 + t4 * t8, t6 * t4, 0, 0, 0, 0; t2 * t8 + t4 * t5, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:49
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (37->9), mult. (41->14), div. (0->0), fcn. (41->6), ass. (0->12)
	t6 = qJ(2) + qJ(3);
	t3 = sin(t6);
	t4 = cos(t6);
	t14 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t17 = t14 + cos(qJ(2)) * pkin(2);
	t15 = r_i_i_C(3) + pkin(8) + pkin(7);
	t13 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t12 = pkin(1) + t17;
	t11 = -sin(qJ(2)) * pkin(2) + t13;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t12 * t8 + t15 * t9, t11 * t9, t13 * t9, 0, 0, 0; t12 * t9 + t15 * t8, t11 * t8, t13 * t8, 0, 0, 0; 0, t17, t14, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:49
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (91->19), mult. (103->23), div. (0->0), fcn. (110->8), ass. (0->20)
	t40 = r_i_i_C(3) + qJ(4);
	t15 = qJ(2) + qJ(3);
	t12 = sin(t15);
	t13 = cos(t15);
	t17 = cos(pkin(11));
	t27 = -r_i_i_C(1) * t17 - pkin(3);
	t16 = sin(pkin(11));
	t34 = r_i_i_C(2) * t16;
	t23 = t40 * t12 + (-t27 - t34) * t13;
	t39 = cos(qJ(2)) * pkin(2) + t23;
	t38 = t12 * t34 + t40 * t13;
	t36 = pkin(1) + t39;
	t19 = sin(qJ(1));
	t30 = t38 * t19;
	t20 = cos(qJ(1));
	t29 = t38 * t20;
	t26 = t27 * t12;
	t24 = t16 * r_i_i_C(1) + t17 * r_i_i_C(2) + pkin(7) + pkin(8);
	t22 = -sin(qJ(2)) * pkin(2) + t26;
	t1 = [-t36 * t19 + t24 * t20, t22 * t20 + t29, t20 * t26 + t29, t20 * t12, 0, 0; t24 * t19 + t36 * t20, t22 * t19 + t30, t19 * t26 + t30, t19 * t12, 0, 0; 0, t39, t23, -t13, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:49
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (152->35), mult. (135->47), div. (0->0), fcn. (146->10), ass. (0->29)
	t11 = cos(pkin(11)) * pkin(4) + pkin(3);
	t20 = qJ(2) + qJ(3);
	t16 = sin(t20);
	t17 = cos(t20);
	t22 = -pkin(9) - qJ(4);
	t41 = t17 * t11 + (r_i_i_C(3) - t22) * t16;
	t18 = cos(qJ(2)) * pkin(2);
	t40 = pkin(1) + t18 + t41;
	t24 = sin(qJ(1));
	t19 = pkin(11) + qJ(5);
	t14 = sin(t19);
	t37 = r_i_i_C(2) * t14;
	t32 = t16 * t37;
	t34 = t17 * t24;
	t39 = r_i_i_C(3) * t34 + t24 * t32;
	t15 = cos(t19);
	t38 = r_i_i_C(1) * t15;
	t25 = cos(qJ(1));
	t33 = t17 * t25;
	t35 = r_i_i_C(3) * t33 + t25 * t32;
	t31 = sin(pkin(11)) * pkin(4) + pkin(8) + pkin(7);
	t29 = -t17 * t22 + (-t11 - t38) * t16;
	t28 = (-t37 + t38) * t17 + t41;
	t27 = -sin(qJ(2)) * pkin(2) + t29;
	t4 = t14 * t24 + t15 * t33;
	t3 = -t14 * t33 + t15 * t24;
	t2 = t14 * t25 - t15 * t34;
	t1 = t14 * t34 + t15 * t25;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t40 * t24 + t31 * t25, t27 * t25 + t35, t29 * t25 + t35, t25 * t16, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t31 * t24 + t40 * t25, t27 * t24 + t39, t29 * t24 + t39, t24 * t16, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, t18 + t28, t28, -t17, (-r_i_i_C(1) * t14 - r_i_i_C(2) * t15) * t16, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:58:49
	% EndTime: 2019-10-10 11:58:50
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (256->43), mult. (177->56), div. (0->0), fcn. (192->12), ass. (0->40)
	t28 = pkin(11) + qJ(5);
	t22 = cos(t28);
	t12 = pkin(5) * t22 + cos(pkin(11)) * pkin(4) + pkin(3);
	t29 = qJ(2) + qJ(3);
	t24 = sin(t29);
	t25 = cos(t29);
	t27 = -pkin(10) - pkin(9) - qJ(4);
	t56 = t25 * t12 + (r_i_i_C(3) - t27) * t24;
	t26 = cos(qJ(2)) * pkin(2);
	t55 = pkin(1) + t26 + t56;
	t23 = qJ(6) + t28;
	t18 = cos(t23);
	t32 = cos(qJ(1));
	t43 = t32 * t18;
	t17 = sin(t23);
	t31 = sin(qJ(1));
	t46 = t31 * t17;
	t5 = t25 * t46 + t43;
	t44 = t32 * t17;
	t45 = t31 * t18;
	t6 = -t25 * t45 + t44;
	t54 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t7 = -t25 * t44 + t45;
	t8 = t25 * t43 + t46;
	t53 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t21 = sin(t28);
	t52 = pkin(5) * t21;
	t51 = r_i_i_C(1) * t18;
	t50 = r_i_i_C(2) * t17;
	t48 = t25 * t31;
	t47 = t25 * t32;
	t39 = t24 * t50;
	t42 = r_i_i_C(3) * t48 + t31 * t39;
	t41 = r_i_i_C(3) * t47 + t32 * t39;
	t40 = t52 + sin(pkin(11)) * pkin(4) + pkin(8) + pkin(7);
	t37 = -r_i_i_C(1) * t17 - r_i_i_C(2) * t18;
	t36 = -t25 * t27 + (-t12 - t51) * t24;
	t35 = (-t50 + t51) * t25 + t56;
	t34 = -sin(qJ(2)) * pkin(2) + t36;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t55 * t31 + t40 * t32, t32 * t34 + t41, t32 * t36 + t41, t32 * t24, (-t21 * t47 + t22 * t31) * pkin(5) + t53, t53; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t40 * t31 + t55 * t32, t31 * t34 + t42, t31 * t36 + t42, t31 * t24, (-t21 * t48 - t22 * t32) * pkin(5) + t54, t54; 0, t26 + t35, t35, -t25, (t37 - t52) * t24, t37 * t24;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end