% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRRR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_jacobia_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(2));
	t1 = sin(qJ(2));
	t3 = [0, -r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; 0, r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t4 = qJ(2) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t6 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t5 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t1 = [0, -sin(qJ(2)) * pkin(2) + t5, t5, 0, 0; 0, cos(qJ(2)) * pkin(2) + t6, t6, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (40->8), mult. (18->8), div. (0->0), fcn. (18->6), ass. (0->9)
	t7 = qJ(2) + qJ(3);
	t6 = qJ(4) + t7;
	t3 = sin(t6);
	t4 = cos(t6);
	t11 = r_i_i_C(1) * t4 - r_i_i_C(2) * t3;
	t10 = t11 + pkin(3) * cos(t7);
	t9 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t8 = -pkin(3) * sin(t7) + t9;
	t1 = [0, -sin(qJ(2)) * pkin(2) + t8, t8, t9, 0; 0, cos(qJ(2)) * pkin(2) + t10, t10, t11, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:05:09
	% EndTime: 2019-12-05 17:05:09
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (87->11), mult. (52->14), div. (0->0), fcn. (52->8), ass. (0->14)
	t23 = pkin(6) + r_i_i_C(3);
	t13 = sin(qJ(5));
	t14 = cos(qJ(5));
	t22 = t14 * r_i_i_C(1) - t13 * r_i_i_C(2);
	t12 = qJ(2) + qJ(3);
	t19 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
	t11 = qJ(4) + t12;
	t8 = sin(t11);
	t9 = cos(t11);
	t18 = -t22 * t8 + t23 * t9;
	t17 = t22 * t9 + t23 * t8;
	t16 = t17 + pkin(3) * cos(t12);
	t15 = -pkin(3) * sin(t12) + t18;
	t1 = [0, -sin(qJ(2)) * pkin(2) + t15, t15, t18, t19 * t9; 0, cos(qJ(2)) * pkin(2) + t16, t16, t17, t19 * t8; 1, 0, 0, 0, t22;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end