% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRP3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRP3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:24:59
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t4 = qJ(1) + qJ(2);
	t2 = sin(t4);
	t3 = cos(t4);
	t6 = r_i_i_C(1) * t3 - r_i_i_C(2) * t2;
	t5 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t1 = [-sin(qJ(1)) * pkin(1) + t5, t5, 0, 0, 0; cos(qJ(1)) * pkin(1) + t6, t6, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (40->8), mult. (18->8), div. (0->0), fcn. (18->6), ass. (0->9)
	t7 = qJ(1) + qJ(2);
	t6 = qJ(3) + t7;
	t3 = sin(t6);
	t4 = cos(t6);
	t11 = r_i_i_C(1) * t4 - r_i_i_C(2) * t3;
	t10 = t11 + pkin(2) * cos(t7);
	t9 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t8 = -pkin(2) * sin(t7) + t9;
	t1 = [-sin(qJ(1)) * pkin(1) + t8, t8, t9, 0, 0; cos(qJ(1)) * pkin(1) + t10, t10, t11, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (105->12), mult. (58->14), div. (0->0), fcn. (58->8), ass. (0->15)
	t25 = r_i_i_C(3) + pkin(8);
	t14 = sin(qJ(4));
	t15 = cos(qJ(4));
	t24 = t15 * r_i_i_C(1) - t14 * r_i_i_C(2);
	t23 = -pkin(3) - t24;
	t13 = qJ(1) + qJ(2);
	t20 = -r_i_i_C(1) * t14 - r_i_i_C(2) * t15;
	t12 = qJ(3) + t13;
	t10 = cos(t12);
	t9 = sin(t12);
	t19 = -t23 * t10 + t25 * t9;
	t18 = t19 + pkin(2) * cos(t13);
	t17 = t25 * t10 + t23 * t9;
	t16 = -pkin(2) * sin(t13) + t17;
	t1 = [-sin(qJ(1)) * pkin(1) + t16, t16, t17, t20 * t10, 0; cos(qJ(1)) * pkin(1) + t18, t18, t19, t20 * t9, 0; 0, 0, 0, t24, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:25:00
	% EndTime: 2019-12-29 20:25:00
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (160->17), mult. (94->17), div. (0->0), fcn. (97->8), ass. (0->18)
	t17 = cos(qJ(4));
	t26 = pkin(4) + r_i_i_C(1);
	t29 = t26 * t17;
	t28 = r_i_i_C(2) + pkin(8);
	t23 = r_i_i_C(3) + qJ(5);
	t16 = sin(qJ(4));
	t27 = t23 * t16 + t29;
	t15 = qJ(1) + qJ(2);
	t14 = qJ(3) + t15;
	t12 = cos(t14);
	t25 = t12 * t16;
	t11 = sin(t14);
	t22 = t28 * t11 + t23 * t25 + (pkin(3) + t29) * t12;
	t21 = pkin(2) * cos(t15) + t22;
	t20 = -t26 * t16 + t23 * t17;
	t19 = (-pkin(3) - t27) * t11 + t28 * t12;
	t18 = -pkin(2) * sin(t15) + t19;
	t1 = [-sin(qJ(1)) * pkin(1) + t18, t18, t19, t20 * t12, t25; cos(qJ(1)) * pkin(1) + t21, t21, t22, t20 * t11, t11 * t16; 0, 0, 0, t27, -t17;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end