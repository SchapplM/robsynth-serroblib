% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RRPR4
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 13:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RRPR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRPR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:33:49
	% EndTime: 2019-12-29 13:33:50
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:33:49
	% EndTime: 2019-12-29 13:33:49
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
	% StartTime: 2019-12-29 13:33:49
	% EndTime: 2019-12-29 13:33:49
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t4 = qJ(1) + qJ(2);
	t2 = sin(t4);
	t3 = cos(t4);
	t6 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t5 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t1 = [-sin(qJ(1)) * pkin(1) + t5, t5, 0, 0; cos(qJ(1)) * pkin(1) + t6, t6, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:33:49
	% EndTime: 2019-12-29 13:33:49
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (41->9), mult. (30->8), div. (0->0), fcn. (32->6), ass. (0->8)
	t18 = qJ(3) + r_i_i_C(3);
	t17 = r_i_i_C(2) * sin(pkin(7)) - pkin(2) - r_i_i_C(1) * cos(pkin(7));
	t10 = qJ(1) + qJ(2);
	t8 = sin(t10);
	t9 = cos(t10);
	t14 = -t17 * t9 + t18 * t8;
	t13 = t17 * t8 + t18 * t9;
	t1 = [-sin(qJ(1)) * pkin(1) + t13, t13, t8, 0; cos(qJ(1)) * pkin(1) + t14, t14, -t9, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 13:33:44
	% EndTime: 2019-12-29 13:33:44
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (70->12), mult. (44->13), div. (0->0), fcn. (46->7), ass. (0->13)
	t21 = pkin(6) + qJ(3) + r_i_i_C(3);
	t11 = pkin(7) + qJ(4);
	t7 = sin(t11);
	t8 = cos(t11);
	t20 = t8 * r_i_i_C(1) - t7 * r_i_i_C(2);
	t19 = -cos(pkin(7)) * pkin(3) - pkin(2) - t20;
	t16 = -r_i_i_C(1) * t7 - r_i_i_C(2) * t8;
	t12 = qJ(1) + qJ(2);
	t10 = cos(t12);
	t9 = sin(t12);
	t15 = -t19 * t10 + t21 * t9;
	t14 = t21 * t10 + t19 * t9;
	t1 = [-sin(qJ(1)) * pkin(1) + t14, t14, t9, t16 * t10; cos(qJ(1)) * pkin(1) + t15, t15, -t10, t16 * t9; 0, 0, 0, t20;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,4);
end