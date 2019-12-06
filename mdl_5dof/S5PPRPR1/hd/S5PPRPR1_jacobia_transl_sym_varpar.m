% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRPR1
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PPRPR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRPR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(pkin(7)), 0, 0, 0; 0, -cos(pkin(7)), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (10->4), mult. (10->6), div. (0->0), fcn. (12->4), ass. (0->7)
	t3 = pkin(8) + qJ(3);
	t1 = sin(t3);
	t2 = cos(t3);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t2;
	t5 = cos(pkin(7));
	t4 = sin(pkin(7));
	t7 = [0, t4, t6 * t5, 0, 0; 0, -t5, t6 * t4, 0, 0; 1, 0, t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (32->8), mult. (33->10), div. (0->0), fcn. (38->6), ass. (0->9)
	t10 = r_i_i_C(3) + qJ(4);
	t9 = r_i_i_C(1) * cos(pkin(9)) - r_i_i_C(2) * sin(pkin(9)) + pkin(3);
	t3 = pkin(8) + qJ(3);
	t1 = sin(t3);
	t2 = cos(t3);
	t8 = -t9 * t1 + t10 * t2;
	t7 = cos(pkin(7));
	t5 = sin(pkin(7));
	t4 = [0, t5, t8 * t7, t7 * t1, 0; 0, -t7, t8 * t5, t5 * t1, 0; 1, 0, t10 * t1 + t9 * t2, -t2, 0;];
	Ja_transl = t4;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:01:48
	% EndTime: 2019-12-05 15:01:48
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (67->17), mult. (56->26), div. (0->0), fcn. (65->7), ass. (0->16)
	t6 = pkin(9) + qJ(5);
	t2 = sin(t6);
	t8 = sin(pkin(7));
	t17 = t8 * t2;
	t4 = cos(t6);
	t16 = t8 * t4;
	t9 = cos(pkin(7));
	t15 = t9 * t2;
	t14 = t9 * t4;
	t13 = r_i_i_C(3) + pkin(6) + qJ(4);
	t12 = t4 * r_i_i_C(1) - t2 * r_i_i_C(2) + cos(pkin(9)) * pkin(4) + pkin(3);
	t7 = pkin(8) + qJ(3);
	t3 = sin(t7);
	t5 = cos(t7);
	t11 = -t12 * t3 + t13 * t5;
	t1 = [0, t8, t11 * t9, t9 * t3, (-t5 * t15 + t16) * r_i_i_C(1) + (-t5 * t14 - t17) * r_i_i_C(2); 0, -t9, t11 * t8, t8 * t3, (-t5 * t17 - t14) * r_i_i_C(1) + (-t5 * t16 + t15) * r_i_i_C(2); 1, 0, t12 * t5 + t13 * t3, -t5, (-r_i_i_C(1) * t2 - r_i_i_C(2) * t4) * t3;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end