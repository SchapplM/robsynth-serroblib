% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RPRR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RPRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:03:01
	% EndTime: 2019-12-29 13:03:01
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:03:01
	% EndTime: 2019-12-29 13:03:01
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:03:01
	% EndTime: 2019-12-29 13:03:01
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(1)) * pkin(1), 0, 0, 0; 0, 1, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:03:01
	% EndTime: 2019-12-29 13:03:01
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (26->8), mult. (12->8), div. (0->0), fcn. (12->6), ass. (0->7)
	t5 = qJ(1) + pkin(7);
	t4 = qJ(3) + t5;
	t2 = sin(t4);
	t3 = cos(t4);
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t2;
	t6 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t1 = [-pkin(2) * sin(t5) - sin(qJ(1)) * pkin(1) + t6, 0, t6, 0; pkin(2) * cos(t5) + cos(qJ(1)) * pkin(1) + t7, 0, t7, 0; 0, 1, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:03:01
	% EndTime: 2019-12-29 13:03:01
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (73->12), mult. (42->14), div. (0->0), fcn. (42->8), ass. (0->13)
	t21 = r_i_i_C(3) + pkin(6);
	t12 = sin(qJ(4));
	t13 = cos(qJ(4));
	t20 = t13 * r_i_i_C(1) - t12 * r_i_i_C(2);
	t19 = -pkin(3) - t20;
	t11 = qJ(1) + pkin(7);
	t16 = -r_i_i_C(1) * t12 - r_i_i_C(2) * t13;
	t10 = qJ(3) + t11;
	t8 = sin(t10);
	t9 = cos(t10);
	t15 = -t19 * t9 + t21 * t8;
	t14 = t19 * t8 + t21 * t9;
	t1 = [-pkin(2) * sin(t11) - sin(qJ(1)) * pkin(1) + t14, 0, t14, t16 * t9; pkin(2) * cos(t11) + cos(qJ(1)) * pkin(1) + t15, 0, t15, t16 * t8; 0, 1, 0, t20;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,4);
end