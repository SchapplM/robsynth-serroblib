% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPPRR1
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:17
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PPPRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPPRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:15
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:15
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(pkin(7)), 0, 0, 0; 0, -cos(pkin(7)), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:15
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (2->2), div. (0->0), fcn. (7->4), ass. (0->4)
	t3 = cos(pkin(7));
	t2 = sin(pkin(7));
	t1 = sin(pkin(8));
	t4 = [0, t2, t3 * t1, 0, 0; 0, -t3, t2 * t1, 0, 0; 1, 0, -cos(pkin(8)), 0, 0;];
	Ja_transl = t4;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (19->10), mult. (22->18), div. (0->0), fcn. (31->6), ass. (0->11)
	t5 = sin(pkin(7));
	t6 = cos(pkin(8));
	t10 = t5 * t6;
	t3 = pkin(9) + qJ(4);
	t1 = sin(t3);
	t7 = cos(pkin(7));
	t9 = t7 * t1;
	t2 = cos(t3);
	t8 = t7 * t2;
	t4 = sin(pkin(8));
	t11 = [0, t5, t7 * t4, (t5 * t2 - t6 * t9) * r_i_i_C(1) + (-t5 * t1 - t6 * t8) * r_i_i_C(2), 0; 0, -t7, t5 * t4, (-t1 * t10 - t8) * r_i_i_C(1) + (-t2 * t10 + t9) * r_i_i_C(2), 0; 1, 0, -t6, (-r_i_i_C(1) * t1 - r_i_i_C(2) * t2) * t4, 0;];
	Ja_transl = t11;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:15
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (72->22), mult. (94->41), div. (0->0), fcn. (119->8), ass. (0->18)
	t19 = pkin(6) + r_i_i_C(3);
	t10 = cos(pkin(8));
	t9 = sin(pkin(7));
	t18 = t10 * t9;
	t12 = sin(qJ(5));
	t8 = sin(pkin(8));
	t17 = t12 * t8;
	t13 = cos(qJ(5));
	t16 = t13 * t8;
	t11 = cos(pkin(7));
	t15 = t10 * t11;
	t14 = r_i_i_C(1) * t13 - r_i_i_C(2) * t12 + pkin(4);
	t7 = pkin(9) + qJ(4);
	t6 = cos(t7);
	t5 = sin(t7);
	t4 = t6 * t15 + t5 * t9;
	t2 = -t11 * t5 + t6 * t18;
	t1 = [0, t9, t11 * t8, t19 * t4 + t14 * (-t5 * t15 + t6 * t9), (t11 * t16 - t4 * t12) * r_i_i_C(1) + (-t11 * t17 - t13 * t4) * r_i_i_C(2); 0, -t11, t9 * t8, t19 * t2 + t14 * (-t11 * t6 - t5 * t18), (-t2 * t12 + t9 * t16) * r_i_i_C(1) + (-t13 * t2 - t9 * t17) * r_i_i_C(2); 1, 0, -t10, (-t14 * t5 + t19 * t6) * t8, (-t10 * t13 - t6 * t17) * r_i_i_C(1) + (t10 * t12 - t6 * t16) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end