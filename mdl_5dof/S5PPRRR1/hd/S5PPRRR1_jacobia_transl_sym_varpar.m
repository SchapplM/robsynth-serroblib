% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRRR1
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
% Datum: 2019-10-24 10:20
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PPRRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:28
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:28
	% EndTime: 2019-10-24 10:20:28
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:28
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(pkin(8)), 0, 0, 0; 0, -cos(pkin(8)), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:29
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (10->4), mult. (10->6), div. (0->0), fcn. (12->4), ass. (0->7)
	t3 = pkin(9) + qJ(3);
	t1 = sin(t3);
	t2 = cos(t3);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t2;
	t5 = cos(pkin(8));
	t4 = sin(pkin(8));
	t7 = [0, t4, t6 * t5, 0, 0; 0, -t5, t6 * t4, 0, 0; 1, 0, t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:28
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (37->7), mult. (25->10), div. (0->0), fcn. (27->6), ass. (0->10)
	t6 = pkin(9) + qJ(3);
	t5 = qJ(4) + t6;
	t2 = sin(t5);
	t3 = cos(t5);
	t11 = r_i_i_C(1) * t3 - r_i_i_C(2) * t2;
	t10 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t9 = -pkin(3) * sin(t6) + t10;
	t8 = cos(pkin(8));
	t7 = sin(pkin(8));
	t1 = [0, t7, t9 * t8, t10 * t8, 0; 0, -t8, t9 * t7, t10 * t7, 0; 1, 0, pkin(3) * cos(t6) + t11, t11, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:20:28
	% EndTime: 2019-10-24 10:20:29
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (110->21), mult. (87->30), div. (0->0), fcn. (93->8), ass. (0->22)
	t15 = pkin(9) + qJ(3);
	t14 = qJ(4) + t15;
	t11 = sin(t14);
	t12 = cos(t14);
	t18 = sin(qJ(5));
	t32 = r_i_i_C(2) * t18;
	t36 = pkin(7) + r_i_i_C(3);
	t37 = t11 * t32 + t12 * t36;
	t19 = cos(qJ(5));
	t35 = -r_i_i_C(1) * t19 - pkin(4);
	t16 = sin(pkin(8));
	t29 = t16 * t18;
	t28 = t16 * t19;
	t17 = cos(pkin(8));
	t27 = t17 * t18;
	t26 = t17 * t19;
	t25 = t37 * t16;
	t24 = t37 * t17;
	t22 = t35 * t11;
	t21 = t36 * t11 + (-t32 - t35) * t12;
	t20 = -pkin(3) * sin(t15) + t22;
	t1 = [0, t16, t20 * t17 + t24, t17 * t22 + t24, (-t12 * t27 + t28) * r_i_i_C(1) + (-t12 * t26 - t29) * r_i_i_C(2); 0, -t17, t20 * t16 + t25, t16 * t22 + t25, (-t12 * t29 - t26) * r_i_i_C(1) + (-t12 * t28 + t27) * r_i_i_C(2); 1, 0, pkin(3) * cos(t15) + t21, t21, (-r_i_i_C(1) * t18 - r_i_i_C(2) * t19) * t11;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end