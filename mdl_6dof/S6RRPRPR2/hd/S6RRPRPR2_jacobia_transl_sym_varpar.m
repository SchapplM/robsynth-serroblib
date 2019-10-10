% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:06
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:08
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
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:08
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
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:08
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (29->10), mult. (31->12), div. (0->0), fcn. (33->6), ass. (0->10)
	t5 = qJ(2) + pkin(10);
	t2 = sin(t5);
	t3 = cos(t5);
	t14 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2) + cos(qJ(2)) * pkin(2);
	t13 = r_i_i_C(3) + qJ(3) + pkin(7);
	t11 = pkin(1) + t14;
	t10 = -sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t11 * t8 + t13 * t9, t10 * t9, t8, 0, 0, 0; t11 * t9 + t13 * t8, t10 * t8, -t9, 0, 0, 0; 0, t14, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:08
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (66->13), mult. (46->16), div. (0->0), fcn. (48->8), ass. (0->13)
	t10 = qJ(2) + pkin(10);
	t7 = qJ(4) + t10;
	t5 = sin(t7);
	t6 = cos(t7);
	t16 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t20 = t16 + pkin(3) * cos(t10) + cos(qJ(2)) * pkin(2);
	t18 = r_i_i_C(3) + pkin(8) + qJ(3) + pkin(7);
	t15 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t14 = pkin(1) + t20;
	t13 = t15 - pkin(3) * sin(t10) - sin(qJ(2)) * pkin(2);
	t12 = cos(qJ(1));
	t11 = sin(qJ(1));
	t1 = [-t14 * t11 + t18 * t12, t13 * t12, t11, t15 * t12, 0, 0; t18 * t11 + t14 * t12, t13 * t11, -t12, t15 * t11, 0, 0; 0, t20, 0, t16, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:08
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (121->21), mult. (76->22), div. (0->0), fcn. (81->8), ass. (0->18)
	t36 = r_i_i_C(3) + qJ(5);
	t18 = qJ(2) + pkin(10);
	t15 = qJ(4) + t18;
	t13 = sin(t15);
	t14 = cos(t15);
	t21 = (pkin(4) - r_i_i_C(2)) * t14 + t36 * t13;
	t35 = t21 + pkin(3) * cos(t18) + cos(qJ(2)) * pkin(2);
	t34 = t36 * t14;
	t33 = pkin(1) + t35;
	t30 = r_i_i_C(1) + pkin(8) + qJ(3) + pkin(7);
	t19 = sin(qJ(1));
	t29 = t19 * t13;
	t20 = cos(qJ(1));
	t28 = t20 * t13;
	t24 = r_i_i_C(2) * t29 + t34 * t19;
	t23 = r_i_i_C(2) * t28 + t34 * t20;
	t22 = -pkin(4) * t13 - pkin(3) * sin(t18) - sin(qJ(2)) * pkin(2);
	t1 = [-t33 * t19 + t30 * t20, t22 * t20 + t23, t19, -pkin(4) * t28 + t23, t28, 0; t30 * t19 + t33 * t20, t22 * t19 + t24, -t20, -pkin(4) * t29 + t24, t29, 0; 0, t35, 0, t21, -t14, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:09
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (192->32), mult. (144->42), div. (0->0), fcn. (157->10), ass. (0->28)
	t26 = sin(qJ(6));
	t28 = cos(qJ(6));
	t51 = r_i_i_C(1) * t26 + r_i_i_C(2) * t28;
	t25 = qJ(2) + pkin(10);
	t22 = qJ(4) + t25;
	t20 = sin(t22);
	t21 = cos(t22);
	t38 = pkin(4) + pkin(9) + r_i_i_C(3);
	t50 = t20 * qJ(5) + t38 * t21;
	t48 = (qJ(5) + t51) * t21;
	t40 = pkin(3) * cos(t25) + cos(qJ(2)) * pkin(2);
	t47 = pkin(1) + t40 + t50;
	t44 = pkin(5) + pkin(8) + qJ(3) + pkin(7);
	t27 = sin(qJ(1));
	t43 = t27 * t26;
	t42 = t27 * t28;
	t29 = cos(qJ(1));
	t41 = t29 * t20;
	t37 = t48 * t27;
	t36 = t48 * t29;
	t32 = t38 * t20;
	t31 = -pkin(3) * sin(t25) - sin(qJ(2)) * pkin(2) - t32;
	t30 = t51 * t20 + t50;
	t4 = -t20 * t43 + t28 * t29;
	t3 = t20 * t42 + t26 * t29;
	t2 = t26 * t41 + t42;
	t1 = t28 * t41 - t43;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t47 * t27 + t44 * t29, t31 * t29 + t36, t27, -t29 * t32 + t36, t41, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t44 * t27 + t47 * t29, t31 * t27 + t37, -t29, -t27 * t32 + t37, t27 * t20, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, t30 + t40, 0, t30, -t21, (-r_i_i_C(1) * t28 + r_i_i_C(2) * t26) * t21;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end