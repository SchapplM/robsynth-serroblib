% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:08:01
	% EndTime: 2020-01-03 12:08:01
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:08:01
	% EndTime: 2020-01-03 12:08:01
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:08:01
	% EndTime: 2020-01-03 12:08:01
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t6 = qJ(1) + qJ(2);
	t4 = sin(t6);
	t5 = cos(t6);
	t8 = t4 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t7 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t1 = [0, 0, 0, 0, 0; cos(qJ(1)) * pkin(1) + t7, t7, 0, 0, 0; sin(qJ(1)) * pkin(1) + t8, t8, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:08:01
	% EndTime: 2020-01-03 12:08:01
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (40->8), mult. (18->8), div. (0->0), fcn. (18->6), ass. (0->9)
	t9 = qJ(1) + qJ(2);
	t8 = qJ(3) + t9;
	t6 = sin(t8);
	t7 = cos(t8);
	t13 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7;
	t12 = pkin(2) * sin(t9) + t13;
	t11 = r_i_i_C(1) * t7 - r_i_i_C(2) * t6;
	t10 = t11 + pkin(2) * cos(t9);
	t1 = [0, 0, 0, 0, 0; cos(qJ(1)) * pkin(1) + t10, t10, t11, 0, 0; sin(qJ(1)) * pkin(1) + t12, t12, t13, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:08:01
	% EndTime: 2020-01-03 12:08:01
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (70->11), mult. (24->10), div. (0->0), fcn. (24->8), ass. (0->10)
	t12 = qJ(1) + qJ(2);
	t11 = qJ(3) + t12;
	t10 = pkin(9) + t11;
	t4 = sin(t10);
	t5 = cos(t10);
	t16 = t4 * r_i_i_C(1) + t5 * r_i_i_C(2) + pkin(3) * sin(t11);
	t15 = pkin(2) * sin(t12) + t16;
	t14 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2) + pkin(3) * cos(t11);
	t13 = t14 + pkin(2) * cos(t12);
	t1 = [0, 0, 0, 1, 0; cos(qJ(1)) * pkin(1) + t13, t13, t14, 0, 0; sin(qJ(1)) * pkin(1) + t15, t15, t16, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:08:01
	% EndTime: 2020-01-03 12:08:01
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (157->16), mult. (64->16), div. (0->0), fcn. (64->10), ass. (0->16)
	t27 = pkin(8) + r_i_i_C(3);
	t16 = sin(qJ(5));
	t17 = cos(qJ(5));
	t26 = t17 * r_i_i_C(1) - t16 * r_i_i_C(2);
	t25 = pkin(4) + t26;
	t15 = qJ(1) + qJ(2);
	t14 = qJ(3) + t15;
	t22 = r_i_i_C(1) * t16 + r_i_i_C(2) * t17;
	t13 = pkin(9) + t14;
	t7 = sin(t13);
	t8 = cos(t13);
	t21 = pkin(3) * cos(t14) + t27 * t7 + t25 * t8;
	t20 = pkin(2) * cos(t15) + t21;
	t19 = pkin(3) * sin(t14) - t27 * t8 + t25 * t7;
	t18 = pkin(2) * sin(t15) + t19;
	t1 = [0, 0, 0, 1, t26; cos(qJ(1)) * pkin(1) + t20, t20, t21, 0, -t22 * t7; sin(qJ(1)) * pkin(1) + t18, t18, t19, 0, t22 * t8;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end