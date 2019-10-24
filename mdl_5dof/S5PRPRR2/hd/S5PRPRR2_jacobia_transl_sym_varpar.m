% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:25
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRPRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:57
	% EndTime: 2019-10-24 10:25:57
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:56
	% EndTime: 2019-10-24 10:25:57
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:56
	% EndTime: 2019-10-24 10:25:57
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t3 = sin(qJ(2));
	t4 = cos(qJ(2));
	t5 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t1 = [0, t5 * cos(pkin(8)), 0, 0, 0; 0, t5 * sin(pkin(8)), 0, 0, 0; 1, r_i_i_C(1) * t4 - r_i_i_C(2) * t3, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:56
	% EndTime: 2019-10-24 10:25:57
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (13->6), mult. (15->8), div. (0->0), fcn. (17->6), ass. (0->7)
	t3 = qJ(2) + pkin(9);
	t1 = sin(t3);
	t2 = cos(t3);
	t7 = -sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t1 - r_i_i_C(2) * t2;
	t5 = cos(pkin(8));
	t4 = sin(pkin(8));
	t6 = [0, t7 * t5, t4, 0, 0; 0, t7 * t4, -t5, 0, 0; 1, t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(2)) * pkin(2), 0, 0, 0;];
	Ja_transl = t6;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:56
	% EndTime: 2019-10-24 10:25:57
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (40->9), mult. (28->12), div. (0->0), fcn. (30->8), ass. (0->10)
	t6 = qJ(2) + pkin(9);
	t5 = qJ(4) + t6;
	t3 = sin(t5);
	t4 = cos(t5);
	t11 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t10 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t9 = -pkin(3) * sin(t6) - sin(qJ(2)) * pkin(2) + t10;
	t8 = cos(pkin(8));
	t7 = sin(pkin(8));
	t1 = [0, t9 * t8, t7, t10 * t8, 0; 0, t9 * t7, -t8, t10 * t7, 0; 1, pkin(3) * cos(t6) + cos(qJ(2)) * pkin(2) + t11, 0, t11, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:25:56
	% EndTime: 2019-10-24 10:25:57
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (113->23), mult. (90->32), div. (0->0), fcn. (96->10), ass. (0->22)
	t15 = qJ(2) + pkin(9);
	t14 = qJ(4) + t15;
	t12 = sin(t14);
	t13 = cos(t14);
	t18 = sin(qJ(5));
	t32 = r_i_i_C(2) * t18;
	t36 = pkin(7) + r_i_i_C(3);
	t37 = t12 * t32 + t13 * t36;
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
	t22 = t35 * t12;
	t21 = -pkin(3) * sin(t15) - sin(qJ(2)) * pkin(2) + t22;
	t20 = t36 * t12 + (-t32 - t35) * t13;
	t1 = [0, t21 * t17 + t24, t16, t17 * t22 + t24, (-t13 * t27 + t28) * r_i_i_C(1) + (-t13 * t26 - t29) * r_i_i_C(2); 0, t21 * t16 + t25, -t17, t16 * t22 + t25, (-t13 * t29 - t26) * r_i_i_C(1) + (-t13 * t28 + t27) * r_i_i_C(2); 1, pkin(3) * cos(t15) + cos(qJ(2)) * pkin(2) + t20, 0, t20, (-r_i_i_C(1) * t18 - r_i_i_C(2) * t19) * t12;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end