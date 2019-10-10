% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:36
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:55
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:55
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
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:55
	% DurationCPUTime: 0.09s
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
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:55
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:55
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (72->13), mult. (51->16), div. (0->0), fcn. (53->8), ass. (0->13)
	t11 = qJ(2) + qJ(3);
	t7 = pkin(10) + t11;
	t4 = sin(t7);
	t5 = cos(t7);
	t18 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2) + pkin(3) * cos(t11);
	t22 = t18 + cos(qJ(2)) * pkin(2);
	t14 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5 - pkin(3) * sin(t11);
	t19 = r_i_i_C(3) + qJ(4) + pkin(8) + pkin(7);
	t16 = pkin(1) + t22;
	t15 = -sin(qJ(2)) * pkin(2) + t14;
	t13 = cos(qJ(1));
	t12 = sin(qJ(1));
	t1 = [-t12 * t16 + t13 * t19, t15 * t13, t14 * t13, t12, 0, 0; t12 * t19 + t13 * t16, t15 * t12, t14 * t12, -t13, 0, 0; 0, t22, t18, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:55
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (165->32), mult. (131->41), div. (0->0), fcn. (141->10), ass. (0->30)
	t24 = qJ(2) + qJ(3);
	t20 = pkin(10) + t24;
	t17 = sin(t20);
	t18 = cos(t20);
	t25 = sin(qJ(5));
	t43 = r_i_i_C(2) * t25;
	t51 = pkin(9) + r_i_i_C(3);
	t52 = t17 * t43 + t18 * t51;
	t49 = t51 * t17 + t18 * pkin(4) + pkin(3) * cos(t24);
	t27 = cos(qJ(5));
	t44 = r_i_i_C(1) * t27;
	t30 = (-pkin(4) - t44) * t17 - pkin(3) * sin(t24);
	t22 = cos(qJ(2)) * pkin(2);
	t47 = pkin(1) + t22 + t49;
	t28 = cos(qJ(1));
	t40 = t25 * t28;
	t26 = sin(qJ(1));
	t39 = t26 * t25;
	t38 = t26 * t27;
	t37 = t27 * t28;
	t36 = t52 * t26;
	t34 = t52 * t28;
	t31 = -sin(qJ(2)) * pkin(2) + t30;
	t29 = (-t43 + t44) * t18 + t49;
	t23 = -qJ(4) - pkin(8) - pkin(7);
	t4 = t18 * t37 + t39;
	t3 = -t18 * t40 + t38;
	t2 = -t18 * t38 + t40;
	t1 = t18 * t39 + t37;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t23 * t28 - t26 * t47, t31 * t28 + t34, t30 * t28 + t34, t26, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t26 * t23 + t28 * t47, t31 * t26 + t36, t30 * t26 + t36, -t28, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, t22 + t29, t29, 0, (-r_i_i_C(1) * t25 - r_i_i_C(2) * t27) * t17, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:36:55
	% EndTime: 2019-10-10 11:36:55
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (247->33), mult. (206->41), div. (0->0), fcn. (229->10), ass. (0->31)
	t28 = sin(qJ(5));
	t30 = cos(qJ(5));
	t38 = r_i_i_C(3) + qJ(6);
	t50 = pkin(5) + r_i_i_C(1);
	t57 = t38 * t28 + t50 * t30;
	t56 = pkin(9) + r_i_i_C(2);
	t27 = qJ(2) + qJ(3);
	t23 = pkin(10) + t27;
	t19 = cos(t23);
	t54 = t19 * t56;
	t18 = sin(t23);
	t53 = t56 * t18 + t19 * pkin(4) + pkin(3) * cos(t27);
	t32 = (-pkin(4) - t57) * t18 - pkin(3) * sin(t27);
	t25 = cos(qJ(2)) * pkin(2);
	t51 = pkin(1) + t25 + t53;
	t29 = sin(qJ(1));
	t48 = t29 * t54;
	t42 = t29 * t28;
	t41 = t29 * t30;
	t31 = cos(qJ(1));
	t40 = t31 * t28;
	t39 = t31 * t30;
	t37 = t31 * t54;
	t35 = t57 * t19 + t53;
	t33 = -sin(qJ(2)) * pkin(2) + t32;
	t26 = -qJ(4) - pkin(8) - pkin(7);
	t4 = t19 * t39 + t42;
	t3 = t19 * t40 - t41;
	t2 = t19 * t41 - t40;
	t1 = t19 * t42 + t39;
	t5 = [-t38 * t1 - t50 * t2 - t31 * t26 - t51 * t29, t33 * t31 + t37, t32 * t31 + t37, t29, -t50 * t3 + t38 * t4, t3; -t29 * t26 + t38 * t3 + t51 * t31 + t50 * t4, t33 * t29 + t48, t32 * t29 + t48, -t31, -t50 * t1 + t38 * t2, t1; 0, t25 + t35, t35, 0, (-t50 * t28 + t38 * t30) * t18, t18 * t28;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end