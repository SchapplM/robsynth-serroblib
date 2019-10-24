% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:38
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPPPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:57
	% EndTime: 2019-10-24 10:38:58
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:57
	% EndTime: 2019-10-24 10:38:57
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [0, 0, 0, 0, 0; -r_i_i_C(1) * t2 + r_i_i_C(2) * t1, 0, 0, 0, 0; -r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:57
	% EndTime: 2019-10-24 10:38:57
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, 1, 0, 0, 0; -t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - cos(qJ(1)) * pkin(1), 0, 0, 0, 0; -t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:57
	% EndTime: 2019-10-24 10:38:57
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (22->8), mult. (16->8), div. (0->0), fcn. (18->6), ass. (0->6)
	t7 = r_i_i_C(3) + qJ(3);
	t6 = -r_i_i_C(1) * cos(pkin(8)) + r_i_i_C(2) * sin(pkin(8)) - pkin(2);
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, 1, 0, 0, 0; -cos(qJ(1)) * pkin(1) - t7 * t1 + t6 * t2, 0, t2, 0, 0; -sin(qJ(1)) * pkin(1) + t7 * t2 + t6 * t1, 0, t1, 0, 0;];
	Ja_transl = t4;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:57
	% EndTime: 2019-10-24 10:38:57
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (42->14), mult. (36->14), div. (0->0), fcn. (45->8), ass. (0->10)
	t4 = sin(pkin(9));
	t6 = cos(pkin(9));
	t9 = t4 * r_i_i_C(1) + t6 * r_i_i_C(2) + qJ(3);
	t5 = sin(pkin(8));
	t7 = cos(pkin(8));
	t8 = -pkin(2) + (-r_i_i_C(3) - qJ(4)) * t5 + (-r_i_i_C(1) * t6 + r_i_i_C(2) * t4 - pkin(3)) * t7;
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [0, 1, 0, -t7, 0; -cos(qJ(1)) * pkin(1) - t9 * t1 + t8 * t2, 0, t2, -t1 * t5, 0; -sin(qJ(1)) * pkin(1) + t9 * t2 + t8 * t1, 0, t1, t2 * t5, 0;];
	Ja_transl = t10;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:38:57
	% EndTime: 2019-10-24 10:38:58
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (83->25), mult. (62->33), div. (0->0), fcn. (75->10), ass. (0->17)
	t14 = cos(pkin(8));
	t11 = qJ(1) + pkin(7);
	t7 = sin(t11);
	t19 = t14 * t7;
	t9 = cos(t11);
	t18 = t14 * t9;
	t17 = pkin(4) * sin(pkin(9)) + qJ(3);
	t13 = sin(pkin(8));
	t16 = -t14 * (cos(pkin(9)) * pkin(4) + pkin(3)) - pkin(2) + (-r_i_i_C(3) - pkin(6) - qJ(4)) * t13;
	t10 = pkin(9) + qJ(5);
	t8 = cos(t10);
	t6 = sin(t10);
	t4 = -t8 * t18 - t6 * t7;
	t3 = t6 * t18 - t7 * t8;
	t2 = t8 * t19 - t6 * t9;
	t1 = t6 * t19 + t8 * t9;
	t5 = [0, 1, 0, -t14, (-r_i_i_C(1) * t6 - r_i_i_C(2) * t8) * t13; -cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t17 * t7 + t16 * t9, 0, t9, -t7 * t13, r_i_i_C(1) * t1 + r_i_i_C(2) * t2; -sin(qJ(1)) * pkin(1) - t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t17 * t9 + t16 * t7, 0, t7, t9 * t13, -r_i_i_C(1) * t3 + r_i_i_C(2) * t4;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end