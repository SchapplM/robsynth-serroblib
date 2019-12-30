% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRPR13_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR13_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR13_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:05:28
	% EndTime: 2019-12-29 17:05:28
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:05:28
	% EndTime: 2019-12-29 17:05:28
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
	% StartTime: 2019-12-29 17:05:28
	% EndTime: 2019-12-29 17:05:28
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:05:28
	% EndTime: 2019-12-29 17:05:28
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (26->9), mult. (26->11), div. (0->0), fcn. (28->5), ass. (0->10)
	t11 = r_i_i_C(3) + pkin(6) + qJ(2);
	t4 = pkin(8) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = cos(pkin(8)) * pkin(2) + pkin(1) + t10;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t8 * t6, t6, t9 * t7, 0, 0; t11 * t6 + t8 * t7, -t7, t9 * t6, 0, 0; 0, 0, t10, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:05:28
	% EndTime: 2019-12-29 17:05:28
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (50->12), mult. (46->13), div. (0->0), fcn. (51->5), ass. (0->12)
	t10 = r_i_i_C(3) + qJ(4);
	t12 = pkin(3) - r_i_i_C(2);
	t4 = pkin(8) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t9 = t10 * t2 + t12 * t3;
	t13 = cos(pkin(8)) * pkin(2) + pkin(1) + t9;
	t11 = r_i_i_C(1) + pkin(6) + qJ(2);
	t8 = t10 * t3 - t12 * t2;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t13 * t6, t6, t8 * t7, t7 * t2, 0; t11 * t6 + t13 * t7, -t7, t8 * t6, t6 * t2, 0; 0, 0, t9, -t3, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 17:05:28
	% EndTime: 2019-12-29 17:05:28
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (87->25), mult. (98->35), div. (0->0), fcn. (111->7), ass. (0->22)
	t18 = pkin(3) + pkin(7) + r_i_i_C(3);
	t8 = pkin(8) + qJ(3);
	t7 = cos(t8);
	t16 = t18 * t7;
	t6 = sin(t8);
	t24 = t16 + t6 * qJ(4) + cos(pkin(8)) * pkin(2) + pkin(1);
	t23 = pkin(4) + pkin(6) + qJ(2);
	t10 = sin(qJ(5));
	t13 = cos(qJ(1));
	t22 = t10 * t13;
	t11 = sin(qJ(1));
	t21 = t11 * t10;
	t12 = cos(qJ(5));
	t20 = t11 * t12;
	t19 = t12 * t13;
	t15 = r_i_i_C(1) * t10 + r_i_i_C(2) * t12 + qJ(4);
	t14 = t15 * t7 - t18 * t6;
	t4 = -t6 * t21 + t19;
	t3 = t6 * t20 + t22;
	t2 = t6 * t22 + t20;
	t1 = t6 * t19 - t21;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t24 * t11 + t23 * t13, t11, t14 * t13, t13 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t23 * t11 + t24 * t13, -t13, t14 * t11, t11 * t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, t15 * t6 + t16, -t7, (-r_i_i_C(1) * t12 + r_i_i_C(2) * t10) * t7;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end