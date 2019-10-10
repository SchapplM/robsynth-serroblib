% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_jacobia_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:36
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:36
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
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:36
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(11);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:36
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (27->9), mult. (26->12), div. (0->0), fcn. (26->6), ass. (0->10)
	t9 = pkin(7) + r_i_i_C(3);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t7 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t6 = pkin(2) + t8;
	t3 = qJ(1) + pkin(11);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [-sin(qJ(1)) * pkin(1) + t9 * t2 - t6 * t1, 0, t7 * t2, 0, 0, 0; cos(qJ(1)) * pkin(1) + t9 * t1 + t6 * t2, 0, t7 * t1, 0, 0, 0; 0, 1, t8, 0, 0, 0;];
	Ja_transl = t10;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:36
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (59->12), mult. (43->16), div. (0->0), fcn. (43->8), ass. (0->13)
	t9 = qJ(3) + qJ(4);
	t5 = sin(t9);
	t6 = cos(t9);
	t15 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t18 = t15 + cos(qJ(3)) * pkin(3);
	t16 = r_i_i_C(3) + pkin(8) + pkin(7);
	t14 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t13 = pkin(2) + t18;
	t12 = -sin(qJ(3)) * pkin(3) + t14;
	t8 = qJ(1) + pkin(11);
	t4 = cos(t8);
	t3 = sin(t8);
	t1 = [-sin(qJ(1)) * pkin(1) + t16 * t4 - t13 * t3, 0, t12 * t4, t14 * t4, 0, 0; cos(qJ(1)) * pkin(1) + t16 * t3 + t13 * t4, 0, t12 * t3, t14 * t3, 0, 0; 0, 1, t18, t15, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:36
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (150->31), mult. (123->44), div. (0->0), fcn. (131->10), ass. (0->26)
	t22 = qJ(3) + qJ(4);
	t18 = sin(t22);
	t19 = cos(t22);
	t23 = sin(qJ(5));
	t42 = pkin(9) + r_i_i_C(3);
	t43 = r_i_i_C(2) * t18 * t23 + t19 * t42;
	t40 = t19 * pkin(4) + t42 * t18;
	t20 = cos(qJ(3)) * pkin(3);
	t39 = pkin(2) + t20 + t40;
	t35 = t19 * t23;
	t25 = cos(qJ(5));
	t34 = t19 * t25;
	t21 = qJ(1) + pkin(11);
	t16 = sin(t21);
	t33 = t43 * t16;
	t17 = cos(t21);
	t32 = t43 * t17;
	t29 = (-r_i_i_C(1) * t25 - pkin(4)) * t18;
	t28 = r_i_i_C(1) * t34 - r_i_i_C(2) * t35 + t40;
	t27 = -sin(qJ(3)) * pkin(3) + t29;
	t26 = -pkin(8) - pkin(7);
	t4 = t16 * t23 + t17 * t34;
	t3 = t16 * t25 - t17 * t35;
	t2 = -t16 * t34 + t17 * t23;
	t1 = t16 * t35 + t17 * t25;
	t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t17 * t26 - t39 * t16, 0, t17 * t27 + t32, t17 * t29 + t32, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t16 * t26 + t39 * t17, 0, t16 * t27 + t33, t16 * t29 + t33, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, 1, t20 + t28, t28, (-r_i_i_C(1) * t23 - r_i_i_C(2) * t25) * t18, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:00:36
	% EndTime: 2019-10-10 09:00:36
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (243->40), mult. (167->55), div. (0->0), fcn. (179->12), ass. (0->35)
	t26 = qJ(5) + qJ(6);
	t20 = sin(t26);
	t27 = qJ(3) + qJ(4);
	t21 = sin(t27);
	t23 = cos(t27);
	t52 = r_i_i_C(2) * t20 * t21 + r_i_i_C(3) * t23;
	t30 = cos(qJ(5));
	t16 = t30 * pkin(5) + pkin(4);
	t31 = -pkin(10) - pkin(9);
	t51 = t23 * t16 + (r_i_i_C(3) - t31) * t21;
	t24 = cos(qJ(3)) * pkin(3);
	t50 = pkin(2) + t24 + t51;
	t25 = qJ(1) + pkin(11);
	t18 = sin(t25);
	t19 = cos(t25);
	t22 = cos(t26);
	t43 = t20 * t23;
	t5 = t18 * t43 + t19 * t22;
	t42 = t22 * t23;
	t6 = -t18 * t42 + t19 * t20;
	t49 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t7 = t18 * t22 - t19 * t43;
	t8 = t18 * t20 + t19 * t42;
	t48 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t28 = sin(qJ(5));
	t47 = pkin(5) * t28;
	t44 = t52 * t18;
	t41 = t23 * t28;
	t40 = t52 * t19;
	t38 = pkin(8) + pkin(7) + t47;
	t36 = -r_i_i_C(1) * t20 - r_i_i_C(2) * t22;
	t35 = -t23 * t31 + (-r_i_i_C(1) * t22 - t16) * t21;
	t34 = r_i_i_C(1) * t42 - r_i_i_C(2) * t43 + t51;
	t33 = -sin(qJ(3)) * pkin(3) + t35;
	t1 = [-sin(qJ(1)) * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t38 * t19 - t50 * t18, 0, t19 * t33 + t40, t19 * t35 + t40, (t18 * t30 - t19 * t41) * pkin(5) + t48, t48; cos(qJ(1)) * pkin(1) + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t38 * t18 + t50 * t19, 0, t18 * t33 + t44, t18 * t35 + t44, (-t18 * t41 - t19 * t30) * pkin(5) + t49, t49; 0, 1, t24 + t34, t34, (t36 - t47) * t21, t36 * t21;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end