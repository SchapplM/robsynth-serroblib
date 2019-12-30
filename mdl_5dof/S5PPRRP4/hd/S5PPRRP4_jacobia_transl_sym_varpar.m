% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRRP4
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:23
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PPRRP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRRP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:23:15
	% EndTime: 2019-12-29 15:23:15
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:23:08
	% EndTime: 2019-12-29 15:23:08
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:23:20
	% EndTime: 2019-12-29 15:23:20
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(pkin(7)), 0, 0, 0; 0, -cos(pkin(7)), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:23:20
	% EndTime: 2019-12-29 15:23:20
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (7->5), mult. (12->8), div. (0->0), fcn. (18->4), ass. (0->7)
	t6 = cos(qJ(3));
	t5 = sin(qJ(3));
	t4 = cos(pkin(7));
	t3 = sin(pkin(7));
	t2 = -t3 * t6 + t4 * t5;
	t1 = -t3 * t5 - t4 * t6;
	t7 = [0, t3, -r_i_i_C(1) * t2 + r_i_i_C(2) * t1, 0, 0; 0, -t4, r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:23:15
	% EndTime: 2019-12-29 15:23:15
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (26->9), mult. (52->14), div. (0->0), fcn. (68->6), ass. (0->13)
	t15 = -pkin(6) - r_i_i_C(3);
	t14 = cos(qJ(3));
	t13 = sin(qJ(3));
	t12 = cos(pkin(7));
	t11 = sin(pkin(7));
	t6 = sin(qJ(4));
	t7 = cos(qJ(4));
	t10 = -t7 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t9 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7;
	t8 = pkin(3) - t10;
	t2 = -t11 * t14 + t12 * t13;
	t1 = -t11 * t13 - t12 * t14;
	t3 = [0, t11, t1 * t15 - t8 * t2, t9 * t1, 0; 0, -t12, t8 * t1 + t15 * t2, t9 * t2, 0; 1, 0, 0, t10, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:23:15
	% EndTime: 2019-12-29 15:23:15
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (37->11), mult. (67->14), div. (0->0), fcn. (89->6), ass. (0->14)
	t17 = pkin(4) + r_i_i_C(1);
	t8 = sin(qJ(4));
	t9 = cos(qJ(4));
	t19 = t8 * r_i_i_C(2) - t17 * t9;
	t16 = -r_i_i_C(3) - qJ(5) - pkin(6);
	t15 = cos(qJ(3));
	t14 = sin(qJ(3));
	t13 = cos(pkin(7));
	t12 = sin(pkin(7));
	t11 = pkin(3) - t19;
	t10 = r_i_i_C(2) * t9 + t17 * t8;
	t2 = -t12 * t15 + t13 * t14;
	t1 = -t12 * t14 - t13 * t15;
	t3 = [0, t12, t1 * t16 - t11 * t2, t10 * t1, t2; 0, -t13, t1 * t11 + t16 * t2, t10 * t2, -t1; 1, 0, 0, t19, 0;];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,5);
end