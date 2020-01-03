% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PRRR6
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4PRRR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:35:18
	% EndTime: 2019-12-31 16:35:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:35:18
	% EndTime: 2019-12-31 16:35:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:35:18
	% EndTime: 2019-12-31 16:35:18
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t3 = sin(qJ(2));
	t4 = cos(qJ(2));
	t5 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t1 = [0, t5 * cos(pkin(7)), 0, 0; 0, t5 * sin(pkin(7)), 0, 0; 1, t4 * r_i_i_C(1) - t3 * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:35:18
	% EndTime: 2019-12-31 16:35:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (19->12), mult. (51->25), div. (0->0), fcn. (55->6), ass. (0->12)
	t3 = sin(qJ(3));
	t6 = cos(qJ(2));
	t11 = t3 * t6;
	t5 = cos(qJ(3));
	t10 = t5 * t6;
	t9 = pkin(5) + r_i_i_C(3);
	t8 = t5 * r_i_i_C(1) - t3 * r_i_i_C(2) + pkin(2);
	t4 = sin(qJ(2));
	t7 = -t8 * t4 + t9 * t6;
	t2 = cos(pkin(7));
	t1 = sin(pkin(7));
	t12 = [0, t7 * t2, (t1 * t5 - t2 * t11) * r_i_i_C(1) + (-t1 * t3 - t2 * t10) * r_i_i_C(2), 0; 0, t7 * t1, (-t1 * t11 - t2 * t5) * r_i_i_C(1) + (-t1 * t10 + t2 * t3) * r_i_i_C(2), 0; 1, t9 * t4 + t8 * t6, (-r_i_i_C(1) * t3 - r_i_i_C(2) * t5) * t4, 0;];
	Ja_transl = t12;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:35:18
	% EndTime: 2019-12-31 16:35:18
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (63->19), mult. (86->35), div. (0->0), fcn. (94->8), ass. (0->19)
	t10 = cos(pkin(7));
	t14 = cos(qJ(2));
	t9 = sin(pkin(7));
	t22 = t14 * t9;
	t8 = qJ(3) + qJ(4);
	t6 = sin(t8);
	t7 = cos(t8);
	t24 = (-t10 * t7 - t6 * t22) * r_i_i_C(1) + (t10 * t6 - t7 * t22) * r_i_i_C(2);
	t20 = t10 * t14;
	t23 = (-t6 * t20 + t7 * t9) * r_i_i_C(1) + (-t7 * t20 - t6 * t9) * r_i_i_C(2);
	t21 = r_i_i_C(3) + pkin(6) + pkin(5);
	t11 = sin(qJ(3));
	t19 = t11 * t14;
	t18 = -r_i_i_C(1) * t6 - r_i_i_C(2) * t7;
	t13 = cos(qJ(3));
	t17 = pkin(3) * t13 + r_i_i_C(1) * t7 - r_i_i_C(2) * t6 + pkin(2);
	t12 = sin(qJ(2));
	t16 = -t17 * t12 + t21 * t14;
	t1 = [0, t16 * t10, (-t10 * t19 + t13 * t9) * pkin(3) + t23, t23; 0, t16 * t9, (-t10 * t13 - t9 * t19) * pkin(3) + t24, t24; 1, t21 * t12 + t17 * t14, (-pkin(3) * t11 + t18) * t12, t18 * t12;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,4);
end