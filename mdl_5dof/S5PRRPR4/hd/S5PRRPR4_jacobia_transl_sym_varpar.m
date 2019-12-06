% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRPR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t3 = sin(qJ(2));
	t4 = cos(qJ(2));
	t5 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t1 = [0, t5 * cos(pkin(8)), 0, 0, 0; 0, t5 * sin(pkin(8)), 0, 0, 0; 1, t4 * r_i_i_C(1) - t3 * r_i_i_C(2), 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (19->12), mult. (51->25), div. (0->0), fcn. (55->6), ass. (0->12)
	t3 = sin(qJ(3));
	t6 = cos(qJ(2));
	t11 = t3 * t6;
	t5 = cos(qJ(3));
	t10 = t5 * t6;
	t9 = pkin(6) + r_i_i_C(3);
	t8 = t5 * r_i_i_C(1) - t3 * r_i_i_C(2) + pkin(2);
	t4 = sin(qJ(2));
	t7 = -t8 * t4 + t9 * t6;
	t2 = cos(pkin(8));
	t1 = sin(pkin(8));
	t12 = [0, t7 * t2, (t1 * t5 - t2 * t11) * r_i_i_C(1) + (-t1 * t3 - t2 * t10) * r_i_i_C(2), 0, 0; 0, t7 * t1, (-t1 * t11 - t2 * t5) * r_i_i_C(1) + (-t1 * t10 + t2 * t3) * r_i_i_C(2), 0, 0; 1, t9 * t4 + t8 * t6, (-r_i_i_C(1) * t3 - r_i_i_C(2) * t5) * t4, 0, 0;];
	Ja_transl = t12;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (47->20), mult. (68->36), div. (0->0), fcn. (75->8), ass. (0->16)
	t17 = r_i_i_C(3) + qJ(4) + pkin(6);
	t11 = cos(qJ(2));
	t5 = sin(pkin(8));
	t16 = t11 * t5;
	t6 = cos(pkin(8));
	t15 = t11 * t6;
	t8 = sin(qJ(3));
	t14 = t11 * t8;
	t10 = cos(qJ(3));
	t4 = qJ(3) + pkin(9);
	t2 = sin(t4);
	t3 = cos(t4);
	t13 = pkin(3) * t10 + r_i_i_C(1) * t3 - r_i_i_C(2) * t2 + pkin(2);
	t9 = sin(qJ(2));
	t12 = t17 * t11 - t13 * t9;
	t1 = [0, t12 * t6, (-t2 * t15 + t5 * t3) * r_i_i_C(1) + (-t3 * t15 - t5 * t2) * r_i_i_C(2) + (t10 * t5 - t6 * t14) * pkin(3), t6 * t9, 0; 0, t12 * t5, (-t2 * t16 - t6 * t3) * r_i_i_C(1) + (-t3 * t16 + t6 * t2) * r_i_i_C(2) + (-t10 * t6 - t5 * t14) * pkin(3), t5 * t9, 0; 1, t13 * t11 + t17 * t9, (-pkin(3) * t8 - r_i_i_C(1) * t2 - r_i_i_C(2) * t3) * t9, -t11, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:24:08
	% EndTime: 2019-12-05 16:24:08
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (109->23), mult. (96->36), div. (0->0), fcn. (107->10), ass. (0->19)
	t14 = qJ(3) + pkin(9);
	t11 = qJ(5) + t14;
	t10 = cos(t11);
	t16 = cos(pkin(8));
	t15 = sin(pkin(8));
	t18 = cos(qJ(2));
	t23 = t15 * t18;
	t9 = sin(t11);
	t26 = (-t16 * t10 - t9 * t23) * r_i_i_C(1) + (-t10 * t23 + t16 * t9) * r_i_i_C(2);
	t22 = t16 * t18;
	t25 = (t15 * t10 - t9 * t22) * r_i_i_C(1) + (-t10 * t22 - t15 * t9) * r_i_i_C(2);
	t24 = r_i_i_C(3) + pkin(7) + qJ(4) + pkin(6);
	t7 = pkin(4) * cos(t14) + cos(qJ(3)) * pkin(3);
	t21 = -r_i_i_C(1) * t9 - r_i_i_C(2) * t10;
	t20 = r_i_i_C(1) * t10 - r_i_i_C(2) * t9 + pkin(2) + t7;
	t17 = sin(qJ(2));
	t19 = -t20 * t17 + t24 * t18;
	t6 = -pkin(4) * sin(t14) - sin(qJ(3)) * pkin(3);
	t1 = [0, t19 * t16, t15 * t7 + t6 * t22 + t25, t16 * t17, t25; 0, t19 * t15, -t16 * t7 + t6 * t23 + t26, t15 * t17, t26; 1, t24 * t17 + t20 * t18, (t21 + t6) * t17, -t18, t21 * t17;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end