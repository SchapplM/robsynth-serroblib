% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRPR3
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
% Datum: 2019-10-24 10:18
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PPRPR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRPR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:50
	% EndTime: 2019-10-24 10:18:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:50
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:50
	% EndTime: 2019-10-24 10:18:50
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(pkin(7)), 0, 0, 0; 0, -cos(pkin(7)), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:50
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->8), mult. (20->15), div. (0->0), fcn. (26->6), ass. (0->10)
	t2 = sin(pkin(7));
	t5 = sin(qJ(3));
	t10 = t2 * t5;
	t6 = cos(qJ(3));
	t9 = t2 * t6;
	t4 = cos(pkin(7));
	t8 = t4 * t5;
	t7 = t4 * t6;
	t3 = cos(pkin(8));
	t1 = [0, t2, (-t3 * t8 + t9) * r_i_i_C(1) + (-t3 * t7 - t10) * r_i_i_C(2), 0, 0; 0, -t4, (-t3 * t10 - t7) * r_i_i_C(1) + (-t3 * t9 + t8) * r_i_i_C(2), 0, 0; 1, 0, (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * sin(pkin(8)), 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:50
	% EndTime: 2019-10-24 10:18:51
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (24->15), mult. (34->26), div. (0->0), fcn. (43->8), ass. (0->14)
	t5 = sin(pkin(7));
	t6 = cos(pkin(8));
	t13 = t5 * t6;
	t8 = sin(qJ(3));
	t12 = t6 * t8;
	t3 = qJ(3) + pkin(9);
	t1 = sin(t3);
	t7 = cos(pkin(7));
	t11 = t7 * t1;
	t2 = cos(t3);
	t10 = t7 * t2;
	t9 = cos(qJ(3));
	t4 = sin(pkin(8));
	t14 = [0, t5, (-t6 * t11 + t5 * t2) * r_i_i_C(1) + (-t5 * t1 - t6 * t10) * r_i_i_C(2) + (-t7 * t12 + t5 * t9) * pkin(3), t7 * t4, 0; 0, -t7, (-t1 * t13 - t10) * r_i_i_C(1) + (-t2 * t13 + t11) * r_i_i_C(2) + (-t5 * t12 - t7 * t9) * pkin(3), t5 * t4, 0; 1, 0, (-pkin(3) * t8 - r_i_i_C(1) * t1 - r_i_i_C(2) * t2) * t4, -t6, 0;];
	Ja_transl = t14;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:50
	% EndTime: 2019-10-24 10:18:50
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (77->27), mult. (106->49), div. (0->0), fcn. (131->10), ass. (0->21)
	t22 = pkin(6) + r_i_i_C(3);
	t10 = cos(pkin(8));
	t9 = sin(pkin(7));
	t21 = t10 * t9;
	t12 = sin(qJ(5));
	t8 = sin(pkin(8));
	t20 = t12 * t8;
	t14 = cos(qJ(5));
	t19 = t14 * t8;
	t11 = cos(pkin(7));
	t18 = t10 * t11;
	t13 = sin(qJ(3));
	t17 = t10 * t13;
	t16 = r_i_i_C(1) * t14 - r_i_i_C(2) * t12 + pkin(4);
	t15 = cos(qJ(3));
	t7 = qJ(3) + pkin(9);
	t6 = cos(t7);
	t5 = sin(t7);
	t4 = t6 * t18 + t5 * t9;
	t2 = -t11 * t5 + t6 * t21;
	t1 = [0, t9, t22 * t4 + (-t11 * t17 + t15 * t9) * pkin(3) + t16 * (-t5 * t18 + t6 * t9), t11 * t8, (t11 * t19 - t12 * t4) * r_i_i_C(1) + (-t11 * t20 - t14 * t4) * r_i_i_C(2); 0, -t11, t22 * t2 + (-t11 * t15 - t9 * t17) * pkin(3) + t16 * (-t11 * t6 - t5 * t21), t9 * t8, (-t12 * t2 + t9 * t19) * r_i_i_C(1) + (-t14 * t2 - t9 * t20) * r_i_i_C(2); 1, 0, (-pkin(3) * t13 - t16 * t5 + t22 * t6) * t8, -t10, (-t10 * t14 - t6 * t20) * r_i_i_C(1) + (t10 * t12 - t6 * t19) * r_i_i_C(2);];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end