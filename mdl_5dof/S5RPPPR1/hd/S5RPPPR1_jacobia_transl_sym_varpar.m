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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
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
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:22
	% EndTime: 2022-01-20 09:13:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:22
	% EndTime: 2022-01-20 09:13:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:22
	% EndTime: 2022-01-20 09:13:22
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(1)) * pkin(1), 0, 0, 0, 0; 0, 1, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:22
	% EndTime: 2022-01-20 09:13:22
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (23->9), mult. (16->8), div. (0->0), fcn. (18->6), ass. (0->6)
	t7 = r_i_i_C(3) + qJ(3);
	t6 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(2);
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-sin(qJ(1)) * pkin(1) + t7 * t2 - t6 * t1, 0, t1, 0, 0; cos(qJ(1)) * pkin(1) + t7 * t1 + t6 * t2, 0, -t2, 0, 0; 0, 1, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:22
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (42->14), mult. (36->14), div. (0->0), fcn. (45->8), ass. (0->10)
	t4 = sin(pkin(9));
	t5 = sin(pkin(8));
	t6 = cos(pkin(9));
	t7 = cos(pkin(8));
	t11 = pkin(2) + (r_i_i_C(3) + qJ(4)) * t5 + (r_i_i_C(1) * t6 - r_i_i_C(2) * t4 + pkin(3)) * t7;
	t8 = t4 * r_i_i_C(1) + t6 * r_i_i_C(2) + qJ(3);
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t9 = [-sin(qJ(1)) * pkin(1) + t8 * t2 - t11 * t1, 0, t1, t2 * t5, 0; cos(qJ(1)) * pkin(1) + t8 * t1 + t11 * t2, 0, -t2, t1 * t5, 0; 0, 1, 0, -t7, 0;];
	Ja_transl = t9;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 09:13:22
	% EndTime: 2022-01-20 09:13:23
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (83->25), mult. (62->33), div. (0->0), fcn. (75->10), ass. (0->17)
	t13 = sin(pkin(8));
	t14 = cos(pkin(8));
	t21 = (r_i_i_C(3) + pkin(6) + qJ(4)) * t13 + t14 * (cos(pkin(9)) * pkin(4) + pkin(3)) + pkin(2);
	t11 = qJ(1) + pkin(7);
	t7 = sin(t11);
	t20 = t14 * t7;
	t9 = cos(t11);
	t19 = t14 * t9;
	t16 = pkin(4) * sin(pkin(9)) + qJ(3);
	t10 = pkin(9) + qJ(5);
	t8 = cos(t10);
	t6 = sin(t10);
	t4 = t8 * t19 + t6 * t7;
	t3 = -t6 * t19 + t7 * t8;
	t2 = -t8 * t20 + t6 * t9;
	t1 = t6 * t20 + t8 * t9;
	t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t9 - t21 * t7, 0, t7, t9 * t13, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t16 * t7 + t21 * t9, 0, -t9, t7 * t13, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 1, 0, -t14, (-r_i_i_C(1) * t6 - r_i_i_C(2) * t8) * t13;];
	Ja_transl = t5;
end