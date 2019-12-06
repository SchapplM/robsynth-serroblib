% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRP1
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRRP1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRP1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:00:21
	% EndTime: 2019-12-05 18:00:22
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:00:21
	% EndTime: 2019-12-05 18:00:22
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:00:22
	% EndTime: 2019-12-05 18:00:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->5), mult. (8->4), div. (0->0), fcn. (10->2), ass. (0->5)
	t4 = pkin(1) - r_i_i_C(2);
	t3 = r_i_i_C(3) + qJ(2);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t4 * t1 + t3 * t2, t1, 0, 0, 0; t3 * t1 + t4 * t2, -t2, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:00:21
	% EndTime: 2019-12-05 18:00:22
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (14->8), mult. (26->10), div. (0->0), fcn. (28->4), ass. (0->9)
	t8 = pkin(1) + pkin(6) + r_i_i_C(3);
	t1 = sin(qJ(3));
	t3 = cos(qJ(3));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = qJ(2) - t6;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t2 * t8 + t4 * t5, t2, t7 * t2, 0, 0; t2 * t5 + t4 * t8, -t4, -t7 * t4, 0, 0; 0, 0, t6, 0, 0;];
	Ja_transl = t9;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:00:21
	% EndTime: 2019-12-05 18:00:22
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (40->14), mult. (45->16), div. (0->0), fcn. (47->6), ass. (0->15)
	t17 = cos(qJ(3)) * pkin(3);
	t5 = qJ(3) + qJ(4);
	t4 = cos(t5);
	t16 = r_i_i_C(1) * t4;
	t3 = sin(t5);
	t15 = r_i_i_C(2) * t3;
	t14 = pkin(1) + r_i_i_C(3) + pkin(7) + pkin(6);
	t13 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t12 = -sin(qJ(3)) * pkin(3) + t13;
	t11 = qJ(2) - t12;
	t9 = cos(qJ(1));
	t7 = sin(qJ(1));
	t2 = t9 * t15;
	t1 = t7 * t16;
	t6 = [t11 * t9 - t14 * t7, t7, t1 + (-t15 + t17) * t7, -t15 * t7 + t1, 0; t11 * t7 + t14 * t9, -t9, t2 + (-t16 - t17) * t9, -t9 * t16 + t2, 0; 0, 0, t12, t13, 0;];
	Ja_transl = t6;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:00:21
	% EndTime: 2019-12-05 18:00:22
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (58->17), mult. (55->18), div. (0->0), fcn. (59->6), ass. (0->17)
	t21 = pkin(4) + r_i_i_C(1);
	t8 = qJ(3) + qJ(4);
	t5 = sin(t8);
	t6 = cos(t8);
	t13 = -t6 * r_i_i_C(2) - t21 * t5;
	t22 = t13 - sin(qJ(3)) * pkin(3);
	t20 = pkin(4) * t6;
	t19 = r_i_i_C(1) * t6;
	t18 = r_i_i_C(2) * t5;
	t14 = pkin(1) + r_i_i_C(3) + qJ(5) + pkin(7) + pkin(6);
	t12 = qJ(2) - t22;
	t11 = cos(qJ(1));
	t10 = sin(qJ(1));
	t4 = t11 * t18;
	t3 = t10 * t19;
	t2 = t20 + cos(qJ(3)) * pkin(3);
	t1 = [-t14 * t10 + t12 * t11, t10, t3 + (t2 - t18) * t10, t3 + (-t18 + t20) * t10, t11; t12 * t10 + t14 * t11, -t11, t4 + (-t2 - t19) * t11, -t11 * t21 * t6 + t4, t10; 0, 0, t22, t13, 0;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end