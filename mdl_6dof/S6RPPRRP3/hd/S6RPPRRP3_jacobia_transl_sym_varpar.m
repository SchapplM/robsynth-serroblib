% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRP3
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:51
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(1)) * pkin(1), 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (19->8), mult. (10->6), div. (0->0), fcn. (12->4), ass. (0->6)
	t5 = pkin(2) - r_i_i_C(2);
	t4 = r_i_i_C(3) + qJ(3);
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t6 = [-sin(qJ(1)) * pkin(1) + t4 * t2 - t5 * t1, 0, t1, 0, 0, 0; cos(qJ(1)) * pkin(1) + t5 * t2 + t4 * t1, 0, -t2, 0, 0, 0; 0, 1, 0, 0, 0, 0;];
	Ja_transl = t6;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (34->11), mult. (28->12), div. (0->0), fcn. (30->6), ass. (0->10)
	t9 = pkin(2) + pkin(7) + r_i_i_C(3);
	t4 = sin(qJ(4));
	t5 = cos(qJ(4));
	t8 = r_i_i_C(1) * t5 - r_i_i_C(2) * t4;
	t7 = -t4 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t6 = qJ(3) - t7;
	t3 = qJ(1) + pkin(9);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [-sin(qJ(1)) * pkin(1) + t6 * t2 - t9 * t1, 0, t1, t8 * t1, 0, 0; cos(qJ(1)) * pkin(1) + t9 * t2 + t6 * t1, 0, -t2, -t8 * t2, 0, 0; 0, 1, 0, t7, 0, 0;];
	Ja_transl = t10;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (80->26), mult. (87->36), div. (0->0), fcn. (97->8), ass. (0->20)
	t11 = cos(qJ(4));
	t15 = pkin(8) + r_i_i_C(3);
	t19 = t15 * t11;
	t10 = cos(qJ(5));
	t8 = sin(qJ(5));
	t13 = r_i_i_C(1) * t10 - r_i_i_C(2) * t8 + pkin(4);
	t9 = sin(qJ(4));
	t18 = t13 * t11 + t15 * t9;
	t17 = t8 * t9;
	t16 = pkin(2) + pkin(7);
	t14 = t10 * t9;
	t12 = t9 * pkin(4) + qJ(3) - t19;
	t7 = qJ(1) + pkin(9);
	t6 = cos(t7);
	t5 = sin(t7);
	t4 = t6 * t14 - t5 * t8;
	t3 = t10 * t5 + t6 * t17;
	t2 = t5 * t14 + t6 * t8;
	t1 = t10 * t6 - t5 * t17;
	t20 = [-sin(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t16 * t5 + t12 * t6, 0, t5, t18 * t5, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; cos(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t6 + t12 * t5, 0, -t6, -t18 * t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 1, 0, -t13 * t9 + t19, (-r_i_i_C(1) * t8 - r_i_i_C(2) * t10) * t11, 0;];
	Ja_transl = t20;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:51:14
	% EndTime: 2019-10-09 23:51:14
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (130->29), mult. (146->37), div. (0->0), fcn. (169->8), ass. (0->22)
	t10 = sin(qJ(4));
	t12 = cos(qJ(4));
	t18 = pkin(8) + r_i_i_C(2);
	t11 = cos(qJ(5));
	t15 = r_i_i_C(3) + qJ(6);
	t19 = pkin(5) + r_i_i_C(1);
	t9 = sin(qJ(5));
	t21 = t19 * t11 + t15 * t9 + pkin(4);
	t24 = t18 * t10 + t21 * t12;
	t22 = t18 * t12;
	t20 = pkin(2) + pkin(7);
	t17 = t10 * t9;
	t16 = t10 * t11;
	t14 = t10 * pkin(4) + qJ(3) - t22;
	t8 = qJ(1) + pkin(9);
	t7 = cos(t8);
	t6 = sin(t8);
	t4 = t7 * t16 - t6 * t9;
	t3 = t11 * t6 + t7 * t17;
	t2 = t6 * t16 + t7 * t9;
	t1 = -t7 * t11 + t6 * t17;
	t5 = [-sin(qJ(1)) * pkin(1) - t20 * t6 + t19 * t4 + t15 * t3 + t14 * t7, 0, t6, t24 * t6, -t19 * t1 + t15 * t2, t1; cos(qJ(1)) * pkin(1) + t20 * t7 + t19 * t2 + t15 * t1 + t14 * t6, 0, -t7, -t24 * t7, -t15 * t4 + t19 * t3, -t3; 0, 1, 0, -t10 * t21 + t22, (t15 * t11 - t19 * t9) * t12, t12 * t9;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end