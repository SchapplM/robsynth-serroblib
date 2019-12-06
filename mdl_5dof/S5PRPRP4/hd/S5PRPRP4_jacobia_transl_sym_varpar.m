% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRP4
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRPRP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:30
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:30
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:31
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (3->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t3 = sin(qJ(2));
	t4 = cos(qJ(2));
	t5 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t1 = [0, t5 * cos(pkin(7)), 0, 0, 0; 0, t5 * sin(pkin(7)), 0, 0, 0; 1, r_i_i_C(1) * t4 - r_i_i_C(2) * t3, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:30
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (13->6), mult. (15->8), div. (0->0), fcn. (17->6), ass. (0->7)
	t3 = qJ(2) + pkin(8);
	t1 = sin(t3);
	t2 = cos(t3);
	t7 = -sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t1 - r_i_i_C(2) * t2;
	t5 = cos(pkin(7));
	t4 = sin(pkin(7));
	t6 = [0, t7 * t5, t4, 0, 0; 0, t7 * t4, -t5, 0, 0; 1, t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(2)) * pkin(2), 0, 0, 0;];
	Ja_transl = t6;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:30
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (44->16), mult. (56->25), div. (0->0), fcn. (62->8), ass. (0->15)
	t4 = sin(pkin(7));
	t6 = sin(qJ(4));
	t15 = t4 * t6;
	t8 = cos(qJ(4));
	t14 = t4 * t8;
	t5 = cos(pkin(7));
	t13 = t5 * t6;
	t12 = t5 * t8;
	t11 = pkin(6) + r_i_i_C(3);
	t10 = t8 * r_i_i_C(1) - t6 * r_i_i_C(2) + pkin(3);
	t3 = qJ(2) + pkin(8);
	t1 = sin(t3);
	t2 = cos(t3);
	t9 = -sin(qJ(2)) * pkin(2) - t10 * t1 + t11 * t2;
	t7 = [0, t9 * t5, t4, (-t2 * t13 + t14) * r_i_i_C(1) + (-t2 * t12 - t15) * r_i_i_C(2), 0; 0, t9 * t4, -t5, (-t2 * t15 - t12) * r_i_i_C(1) + (-t2 * t14 + t13) * r_i_i_C(2), 0; 1, cos(qJ(2)) * pkin(2) + t11 * t1 + t10 * t2, 0, (-r_i_i_C(1) * t6 - r_i_i_C(2) * t8) * t1, 0;];
	Ja_transl = t7;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:36:30
	% EndTime: 2019-12-05 15:36:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (77->18), mult. (99->26), div. (0->0), fcn. (114->8), ass. (0->19)
	t11 = sin(qJ(4));
	t13 = cos(qJ(4));
	t15 = r_i_i_C(3) + qJ(5);
	t21 = pkin(4) + r_i_i_C(1);
	t22 = t15 * t11 + t21 * t13 + pkin(3);
	t20 = pkin(6) + r_i_i_C(2);
	t9 = sin(pkin(7));
	t19 = t9 * t11;
	t18 = t9 * t13;
	t10 = cos(pkin(7));
	t17 = t10 * t11;
	t16 = t10 * t13;
	t8 = qJ(2) + pkin(8);
	t6 = sin(t8);
	t7 = cos(t8);
	t14 = -sin(qJ(2)) * pkin(2) + t20 * t7 - t22 * t6;
	t3 = t7 * t17 - t18;
	t1 = t7 * t19 + t16;
	t2 = [0, t14 * t10, t9, t15 * (t7 * t16 + t19) - t21 * t3, t3; 0, t14 * t9, -t10, t15 * (t7 * t18 - t17) - t21 * t1, t1; 1, cos(qJ(2)) * pkin(2) + t20 * t6 + t22 * t7, 0, (-t21 * t11 + t15 * t13) * t6, t6 * t11;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,5);
end