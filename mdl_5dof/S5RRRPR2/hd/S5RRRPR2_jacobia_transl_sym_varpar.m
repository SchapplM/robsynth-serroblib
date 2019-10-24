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
% Datum: 2019-10-24 10:50
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:50:34
	% EndTime: 2019-10-24 10:50:34
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:34
	% EndTime: 2019-10-24 10:50:34
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [0, 0, 0, 0, 0; -r_i_i_C(1) * t2 + r_i_i_C(2) * t1, 0, 0, 0, 0; -r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:34
	% EndTime: 2019-10-24 10:50:34
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t4 = qJ(1) + qJ(2);
	t2 = sin(t4);
	t3 = cos(t4);
	t6 = -t3 * r_i_i_C(1) + t2 * r_i_i_C(2);
	t5 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t1 = [0, 0, 0, 0, 0; -cos(qJ(1)) * pkin(1) + t6, t6, 0, 0, 0; -sin(qJ(1)) * pkin(1) + t5, t5, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:34
	% EndTime: 2019-10-24 10:50:34
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (40->8), mult. (18->8), div. (0->0), fcn. (18->6), ass. (0->9)
	t7 = qJ(1) + qJ(2);
	t6 = qJ(3) + t7;
	t2 = sin(t6);
	t3 = cos(t6);
	t11 = -t3 * r_i_i_C(1) + t2 * r_i_i_C(2);
	t10 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t9 = -pkin(2) * cos(t7) + t11;
	t8 = -pkin(2) * sin(t7) + t10;
	t1 = [0, 0, 0, 0, 0; -cos(qJ(1)) * pkin(1) + t9, t9, t11, 0, 0; -sin(qJ(1)) * pkin(1) + t8, t8, t10, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:34
	% EndTime: 2019-10-24 10:50:34
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (70->11), mult. (24->10), div. (0->0), fcn. (24->8), ass. (0->10)
	t10 = qJ(1) + qJ(2);
	t9 = qJ(3) + t10;
	t4 = pkin(9) + t9;
	t2 = sin(t4);
	t3 = cos(t4);
	t14 = -pkin(3) * cos(t9) - t3 * r_i_i_C(1) + t2 * r_i_i_C(2);
	t13 = -pkin(3) * sin(t9) - t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t12 = -pkin(2) * cos(t10) + t14;
	t11 = -pkin(2) * sin(t10) + t13;
	t1 = [0, 0, 0, 1, 0; -cos(qJ(1)) * pkin(1) + t12, t12, t14, 0, 0; -sin(qJ(1)) * pkin(1) + t11, t11, t13, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:50:34
	% EndTime: 2019-10-24 10:50:34
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (157->16), mult. (64->16), div. (0->0), fcn. (64->10), ass. (0->16)
	t14 = sin(qJ(5));
	t15 = cos(qJ(5));
	t26 = r_i_i_C(1) * t15 - t14 * r_i_i_C(2);
	t25 = pkin(8) + r_i_i_C(3);
	t24 = -pkin(4) - t26;
	t13 = qJ(1) + qJ(2);
	t12 = qJ(3) + t13;
	t20 = r_i_i_C(1) * t14 + r_i_i_C(2) * t15;
	t7 = pkin(9) + t12;
	t5 = sin(t7);
	t6 = cos(t7);
	t19 = -pkin(3) * sin(t12) + t25 * t6 + t24 * t5;
	t18 = -pkin(3) * cos(t12) - t25 * t5 + t24 * t6;
	t17 = -pkin(2) * sin(t13) + t19;
	t16 = -pkin(2) * cos(t13) + t18;
	t1 = [0, 0, 0, 1, t26; -cos(qJ(1)) * pkin(1) + t16, t16, t18, 0, t20 * t5; -sin(qJ(1)) * pkin(1) + t17, t17, t19, 0, -t20 * t6;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end