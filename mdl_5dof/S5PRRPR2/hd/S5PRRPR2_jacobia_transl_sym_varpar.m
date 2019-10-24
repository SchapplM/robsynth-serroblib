% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPR2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:30
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:15
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:15
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:15
	% EndTime: 2019-10-24 10:30:15
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (6->3), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->4)
	t3 = pkin(8) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, -t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0, 0; 0, t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:15
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (24->6), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->7)
	t5 = pkin(8) + qJ(2);
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
	% StartTime: 2019-10-24 10:30:15
	% EndTime: 2019-10-24 10:30:15
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (65->10), mult. (30->8), div. (0->0), fcn. (32->6), ass. (0->9)
	t19 = qJ(4) + r_i_i_C(3);
	t18 = r_i_i_C(2) * sin(pkin(9)) - pkin(3) - r_i_i_C(1) * cos(pkin(9));
	t11 = pkin(8) + qJ(2);
	t10 = qJ(3) + t11;
	t8 = sin(t10);
	t9 = cos(t10);
	t15 = -t18 * t9 + t19 * t8;
	t14 = t18 * t8 + t19 * t9;
	t1 = [0, -pkin(2) * sin(t11) + t14, t14, t8, 0; 0, pkin(2) * cos(t11) + t15, t15, -t9, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:15
	% EndTime: 2019-10-24 10:30:16
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (136->21), mult. (86->29), div. (0->0), fcn. (100->8), ass. (0->18)
	t19 = sin(pkin(9));
	t20 = cos(pkin(9));
	t31 = pkin(4) * t20 + (pkin(7) + r_i_i_C(3)) * t19 + pkin(3);
	t21 = sin(qJ(5));
	t26 = t20 * t21;
	t22 = cos(qJ(5));
	t25 = t20 * t22;
	t18 = pkin(8) + qJ(2);
	t17 = qJ(3) + t18;
	t15 = sin(t17);
	t16 = cos(t17);
	t7 = t15 * t22 - t16 * t26;
	t8 = t15 * t21 + t16 * t25;
	t24 = t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t15 * qJ(4) + t31 * t16;
	t5 = t15 * t26 + t16 * t22;
	t6 = -t15 * t25 + t16 * t21;
	t23 = t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t16 * qJ(4) - t31 * t15;
	t1 = [0, -pkin(2) * sin(t18) + t23, t23, t15, r_i_i_C(1) * t7 - r_i_i_C(2) * t8; 0, pkin(2) * cos(t18) + t24, t24, -t16, -t5 * r_i_i_C(1) + t6 * r_i_i_C(2); 1, 0, 0, 0, (-r_i_i_C(1) * t21 - r_i_i_C(2) * t22) * t19;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end