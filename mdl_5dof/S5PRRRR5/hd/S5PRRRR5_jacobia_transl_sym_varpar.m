% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR5
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR5_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR5_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:43:07
	% EndTime: 2019-12-29 15:43:07
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:43:12
	% EndTime: 2019-12-29 15:43:12
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:43:07
	% EndTime: 2019-12-29 15:43:07
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (6->3), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->4)
	t3 = pkin(9) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, -t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0, 0; 0, t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:43:07
	% EndTime: 2019-12-29 15:43:07
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (25->7), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->10)
	t9 = pkin(6) + r_i_i_C(3);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t7 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t6 = pkin(2) + t8;
	t3 = pkin(9) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [0, -t6 * t1 + t9 * t2, t7 * t2, 0, 0; 0, t9 * t1 + t6 * t2, t7 * t1, 0, 0; 1, 0, t8, 0, 0;];
	Ja_transl = t10;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:43:12
	% EndTime: 2019-12-29 15:43:12
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (57->10), mult. (41->14), div. (0->0), fcn. (41->6), ass. (0->13)
	t9 = qJ(3) + qJ(4);
	t5 = sin(t9);
	t6 = cos(t9);
	t15 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t18 = t15 + cos(qJ(3)) * pkin(3);
	t16 = r_i_i_C(3) + pkin(7) + pkin(6);
	t14 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t13 = pkin(2) + t18;
	t12 = -sin(qJ(3)) * pkin(3) + t14;
	t8 = pkin(9) + qJ(2);
	t4 = cos(t8);
	t3 = sin(t8);
	t1 = [0, -t13 * t3 + t16 * t4, t12 * t4, t14 * t4, 0; 0, t13 * t4 + t16 * t3, t12 * t3, t14 * t3, 0; 1, 0, t18, t15, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:43:07
	% EndTime: 2019-12-29 15:43:07
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (112->13), mult. (61->18), div. (0->0), fcn. (61->8), ass. (0->16)
	t14 = qJ(3) + qJ(4);
	t10 = qJ(5) + t14;
	t5 = sin(t10);
	t6 = cos(t10);
	t20 = t6 * r_i_i_C(1) - r_i_i_C(2) * t5;
	t19 = t20 + pkin(4) * cos(t14);
	t24 = cos(qJ(3)) * pkin(3) + t19;
	t18 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t15 = t18 - pkin(4) * sin(t14);
	t21 = r_i_i_C(3) + pkin(8) + pkin(7) + pkin(6);
	t17 = pkin(2) + t24;
	t16 = -sin(qJ(3)) * pkin(3) + t15;
	t12 = pkin(9) + qJ(2);
	t8 = cos(t12);
	t7 = sin(t12);
	t1 = [0, -t17 * t7 + t21 * t8, t16 * t8, t15 * t8, t18 * t8; 0, t17 * t8 + t21 * t7, t16 * t7, t15 * t7, t18 * t7; 1, 0, t24, t19, t20;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end