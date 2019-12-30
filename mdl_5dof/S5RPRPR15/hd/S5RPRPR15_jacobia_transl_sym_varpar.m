% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR15
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRPR15_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR15_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR15_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:10:50
	% EndTime: 2019-12-29 17:10:50
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:10:50
	% EndTime: 2019-12-29 17:10:50
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:10:45
	% EndTime: 2019-12-29 17:10:45
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (7->5), mult. (8->4), div. (0->0), fcn. (10->2), ass. (0->5)
	t4 = pkin(1) - r_i_i_C(2);
	t3 = r_i_i_C(3) + qJ(2);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t1 * t4 + t2 * t3, t1, 0, 0, 0; t1 * t3 + t2 * t4, -t2, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:10:45
	% EndTime: 2019-12-29 17:10:45
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (14->8), mult. (26->10), div. (0->0), fcn. (28->4), ass. (0->9)
	t8 = pkin(1) + pkin(6) + r_i_i_C(3);
	t1 = sin(qJ(3));
	t3 = cos(qJ(3));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -t1 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t5 = qJ(2) - t6;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t8 * t2 + t5 * t4, t2, t7 * t2, 0, 0; t5 * t2 + t8 * t4, -t4, -t7 * t4, 0, 0; 0, 0, t6, 0, 0;];
	Ja_transl = t9;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:10:45
	% EndTime: 2019-12-29 17:10:45
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (32->14), mult. (67->16), div. (0->0), fcn. (76->6), ass. (0->13)
	t10 = r_i_i_C(3) + qJ(4);
	t3 = sin(qJ(3));
	t5 = cos(qJ(3));
	t1 = sin(pkin(8));
	t2 = cos(pkin(8));
	t9 = r_i_i_C(1) * t2 - r_i_i_C(2) * t1 + pkin(3);
	t12 = t10 * t5 - t9 * t3;
	t11 = t10 * t3 + t9 * t5;
	t8 = r_i_i_C(1) * t1 + r_i_i_C(2) * t2 + pkin(1) + pkin(6);
	t7 = qJ(2) - t12;
	t6 = cos(qJ(1));
	t4 = sin(qJ(1));
	t13 = [-t8 * t4 + t7 * t6, t4, t11 * t4, -t4 * t5, 0; t7 * t4 + t8 * t6, -t6, -t11 * t6, t6 * t5, 0; 0, 0, t12, t3, 0;];
	Ja_transl = t13;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:10:50
	% EndTime: 2019-12-29 17:10:50
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (75->27), mult. (96->38), div. (0->0), fcn. (109->8), ass. (0->21)
	t13 = cos(qJ(3));
	t20 = r_i_i_C(3) + pkin(7) + qJ(4);
	t22 = t20 * t13;
	t11 = sin(qJ(3));
	t5 = cos(pkin(8)) * pkin(4) + pkin(3);
	t8 = pkin(8) + qJ(5);
	t6 = sin(t8);
	t7 = cos(t8);
	t16 = r_i_i_C(1) * t7 - r_i_i_C(2) * t6 + t5;
	t21 = t20 * t11 + t16 * t13;
	t12 = sin(qJ(1));
	t19 = t11 * t12;
	t14 = cos(qJ(1));
	t18 = t11 * t14;
	t17 = pkin(4) * sin(pkin(8)) + pkin(1) + pkin(6);
	t15 = t11 * t5 + qJ(2) - t22;
	t4 = -t12 * t6 + t7 * t18;
	t3 = t12 * t7 + t6 * t18;
	t2 = t14 * t6 + t7 * t19;
	t1 = t14 * t7 - t6 * t19;
	t9 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t17 * t12 + t15 * t14, t12, t21 * t12, -t12 * t13, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * t12 + t17 * t14, -t14, -t21 * t14, t14 * t13, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, -t16 * t11 + t22, t11, (-r_i_i_C(1) * t6 - r_i_i_C(2) * t7) * t13;];
	Ja_transl = t9;
else
	Ja_transl=NaN(3,5);
end