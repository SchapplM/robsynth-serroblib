% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR3
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
% Datum: 2019-10-24 10:36
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRRR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:13
	% EndTime: 2019-10-24 10:36:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:13
	% EndTime: 2019-10-24 10:36:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:13
	% EndTime: 2019-10-24 10:36:13
	% DurationCPUTime: 0.06s
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
	% StartTime: 2019-10-24 10:36:13
	% EndTime: 2019-10-24 10:36:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (24->6), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->7)
	t5 = pkin(9) + qJ(2);
	t4 = qJ(3) + t5;
	t2 = sin(t4);
	t3 = cos(t4);
	t7 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t6 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t1 = [0, -pkin(2) * sin(t5) + t6, t6, 0, 0; 0, pkin(2) * cos(t5) + t7, t7, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:13
	% EndTime: 2019-10-24 10:36:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (58->9), mult. (18->8), div. (0->0), fcn. (18->6), ass. (0->10)
	t8 = pkin(9) + qJ(2);
	t7 = qJ(3) + t8;
	t6 = qJ(4) + t7;
	t3 = sin(t6);
	t4 = cos(t6);
	t12 = r_i_i_C(1) * t4 - r_i_i_C(2) * t3;
	t11 = t12 + pkin(3) * cos(t7);
	t10 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t9 = -pkin(3) * sin(t7) + t10;
	t1 = [0, -pkin(2) * sin(t8) + t9, t9, t10, 0; 0, pkin(2) * cos(t8) + t11, t11, t12, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:13
	% EndTime: 2019-10-24 10:36:13
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (145->13), mult. (58->14), div. (0->0), fcn. (58->8), ass. (0->16)
	t26 = r_i_i_C(3) + pkin(8);
	t15 = sin(qJ(5));
	t16 = cos(qJ(5));
	t25 = t16 * r_i_i_C(1) - t15 * r_i_i_C(2);
	t24 = -pkin(4) - t25;
	t14 = pkin(9) + qJ(2);
	t13 = qJ(3) + t14;
	t21 = -r_i_i_C(1) * t15 - r_i_i_C(2) * t16;
	t12 = qJ(4) + t13;
	t10 = cos(t12);
	t9 = sin(t12);
	t20 = -t24 * t10 + t26 * t9;
	t19 = t20 + pkin(3) * cos(t13);
	t18 = t26 * t10 + t24 * t9;
	t17 = -pkin(3) * sin(t13) + t18;
	t1 = [0, -pkin(2) * sin(t14) + t17, t17, t18, t21 * t10; 0, pkin(2) * cos(t14) + t19, t19, t20, t21 * t9; 1, 0, 0, 0, t25;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end