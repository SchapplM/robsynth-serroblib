% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPPRP5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRP5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:05:12
	% EndTime: 2019-12-29 16:05:12
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:05:12
	% EndTime: 2019-12-29 16:05:12
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:05:06
	% EndTime: 2019-12-29 16:05:07
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(7)) - r_i_i_C(2) * sin(pkin(7)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:05:12
	% EndTime: 2019-12-29 16:05:12
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (14->9), mult. (24->8), div. (0->0), fcn. (29->4), ass. (0->7)
	t1 = sin(pkin(7));
	t2 = cos(pkin(7));
	t8 = pkin(1) + (pkin(2) + r_i_i_C(1)) * t2 + (r_i_i_C(3) + qJ(3)) * t1;
	t6 = r_i_i_C(2) + qJ(2);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t5 = [-t3 * t8 + t6 * t4, t3, t4 * t1, 0, 0; t6 * t3 + t4 * t8, -t4, t3 * t1, 0, 0; 0, 0, -t2, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:05:12
	% EndTime: 2019-12-29 16:05:12
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (31->17), mult. (68->26), div. (0->0), fcn. (83->6), ass. (0->15)
	t5 = sin(pkin(7));
	t6 = cos(pkin(7));
	t16 = qJ(3) * t5 + pkin(1) + (pkin(2) + pkin(3)) * t6;
	t14 = -pkin(6) - r_i_i_C(3) + qJ(2);
	t7 = sin(qJ(4));
	t9 = cos(qJ(4));
	t12 = t5 * t9 - t6 * t7;
	t11 = t5 * t7 + t6 * t9;
	t10 = cos(qJ(1));
	t8 = sin(qJ(1));
	t4 = t11 * t10;
	t3 = t12 * t10;
	t2 = t11 * t8;
	t1 = t12 * t8;
	t13 = [-t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + t14 * t10 - t16 * t8, t8, t10 * t5, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t16 * t10 + t14 * t8, -t10, t8 * t5, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; 0, 0, -t6, -t11 * r_i_i_C(1) - t12 * r_i_i_C(2), 0;];
	Ja_transl = t13;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:05:12
	% EndTime: 2019-12-29 16:05:12
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (54->20), mult. (124->27), div. (0->0), fcn. (155->6), ass. (0->19)
	t13 = cos(qJ(4));
	t10 = cos(pkin(7));
	t11 = sin(qJ(4));
	t18 = t10 * t11;
	t9 = sin(pkin(7));
	t24 = -t13 * t9 + t18;
	t23 = qJ(3) * t9 + pkin(1) + (pkin(2) + pkin(3)) * t10;
	t21 = pkin(4) + r_i_i_C(1);
	t14 = cos(qJ(1));
	t19 = t14 * t9;
	t17 = r_i_i_C(3) + qJ(5);
	t16 = -pkin(6) - r_i_i_C(2) + qJ(2);
	t5 = t10 * t13 + t11 * t9;
	t12 = sin(qJ(1));
	t4 = t5 * t14;
	t3 = -t13 * t19 + t14 * t18;
	t2 = t5 * t12;
	t1 = t24 * t12;
	t6 = [-t17 * t1 - t23 * t12 + t16 * t14 - t21 * t2, t12, t19, t17 * t4 - t21 * t3, t3; t16 * t12 + t23 * t14 + t17 * t3 + t21 * t4, -t14, t12 * t9, -t21 * t1 + t17 * t2, t1; 0, 0, -t10, -t17 * t24 - t21 * t5, t5;];
	Ja_transl = t6;
else
	Ja_transl=NaN(3,5);
end