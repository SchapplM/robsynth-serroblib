% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:21
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PPRRR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRRR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(pkin(8)), 0, 0, 0; 0, -cos(pkin(8)), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->8), mult. (20->15), div. (0->0), fcn. (26->6), ass. (0->10)
	t2 = sin(pkin(8));
	t5 = sin(qJ(3));
	t10 = t2 * t5;
	t6 = cos(qJ(3));
	t9 = t2 * t6;
	t4 = cos(pkin(8));
	t8 = t4 * t5;
	t7 = t4 * t6;
	t3 = cos(pkin(9));
	t1 = [0, t2, (-t3 * t8 + t9) * r_i_i_C(1) + (-t3 * t7 - t10) * r_i_i_C(2), 0, 0; 0, -t4, (-t3 * t10 - t7) * r_i_i_C(1) + (-t3 * t9 + t8) * r_i_i_C(2), 0, 0; 1, 0, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * sin(pkin(9)), 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:17
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (36->20), mult. (92->38), div. (0->0), fcn. (114->8), ass. (0->20)
	t5 = sin(pkin(9));
	t9 = sin(qJ(4));
	t21 = t5 * t9;
	t20 = pkin(6) + r_i_i_C(3);
	t11 = cos(qJ(4));
	t19 = t11 * t5;
	t12 = cos(qJ(3));
	t18 = t12 * t5;
	t6 = sin(pkin(8));
	t17 = t12 * t6;
	t8 = cos(pkin(8));
	t16 = t12 * t8;
	t10 = sin(qJ(3));
	t15 = t6 * t10;
	t14 = t8 * t10;
	t13 = r_i_i_C(1) * t11 - r_i_i_C(2) * t9 + pkin(3);
	t7 = cos(pkin(9));
	t4 = t7 * t16 + t15;
	t2 = t7 * t17 - t14;
	t1 = [0, t6, t20 * t4 + t13 * (-t7 * t14 + t17), (t8 * t19 - t4 * t9) * r_i_i_C(1) + (-t11 * t4 - t8 * t21) * r_i_i_C(2), 0; 0, -t8, t20 * t2 + t13 * (-t7 * t15 - t16), (t6 * t19 - t2 * t9) * r_i_i_C(1) + (-t11 * t2 - t6 * t21) * r_i_i_C(2), 0; 1, 0, (-t13 * t10 + t20 * t12) * t5, (-t7 * t11 - t9 * t18) * r_i_i_C(1) + (-t11 * t18 + t7 * t9) * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:21:17
	% EndTime: 2019-10-24 10:21:18
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (93->28), mult. (152->49), div. (0->0), fcn. (186->10), ass. (0->27)
	t14 = qJ(4) + qJ(5);
	t12 = sin(t14);
	t13 = cos(t14);
	t15 = sin(pkin(9));
	t16 = sin(pkin(8));
	t32 = t15 * t16;
	t17 = cos(pkin(9));
	t18 = cos(pkin(8));
	t20 = sin(qJ(3));
	t26 = t18 * t20;
	t22 = cos(qJ(3));
	t27 = t16 * t22;
	t8 = t17 * t27 - t26;
	t36 = (-t12 * t8 + t13 * t32) * r_i_i_C(1) + (-t12 * t32 - t13 * t8) * r_i_i_C(2);
	t25 = t18 * t22;
	t28 = t16 * t20;
	t10 = t17 * t25 + t28;
	t31 = t15 * t18;
	t35 = (-t10 * t12 + t13 * t31) * r_i_i_C(1) + (-t10 * t13 - t12 * t31) * r_i_i_C(2);
	t29 = t15 * t22;
	t34 = (-t12 * t29 - t13 * t17) * r_i_i_C(1) + (t12 * t17 - t13 * t29) * r_i_i_C(2);
	t33 = r_i_i_C(3) + pkin(7) + pkin(6);
	t21 = cos(qJ(4));
	t30 = t15 * t21;
	t24 = pkin(4) * t21 + r_i_i_C(1) * t13 - r_i_i_C(2) * t12 + pkin(3);
	t19 = sin(qJ(4));
	t1 = [0, t16, t33 * t10 + t24 * (-t17 * t26 + t27), (-t10 * t19 + t18 * t30) * pkin(4) + t35, t35; 0, -t18, t33 * t8 + t24 * (-t17 * t28 - t25), (t16 * t30 - t19 * t8) * pkin(4) + t36, t36; 1, 0, (-t24 * t20 + t33 * t22) * t15, (-t17 * t21 - t19 * t29) * pkin(4) + t34, t34;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end