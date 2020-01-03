% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPPR5
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRPPR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPPR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPPR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:38:49
	% EndTime: 2019-12-31 17:38:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:38:50
	% EndTime: 2019-12-31 17:38:50
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:38:50
	% EndTime: 2019-12-31 17:38:50
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t3 = sin(qJ(2));
	t4 = cos(qJ(2));
	t5 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t1 = [0, t5 * cos(pkin(7)), 0, 0, 0; 0, t5 * sin(pkin(7)), 0, 0, 0; 1, t4 * r_i_i_C(1) - t3 * r_i_i_C(2), 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:38:50
	% EndTime: 2019-12-31 17:38:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (10->5), mult. (22->8), div. (0->0), fcn. (25->4), ass. (0->8)
	t7 = pkin(2) + r_i_i_C(1);
	t6 = r_i_i_C(3) + qJ(3);
	t3 = sin(qJ(2));
	t4 = cos(qJ(2));
	t5 = -t3 * t7 + t4 * t6;
	t2 = cos(pkin(7));
	t1 = sin(pkin(7));
	t8 = [0, t5 * t2, t2 * t3, 0, 0; 0, t5 * t1, t1 * t3, 0, 0; 1, t3 * t6 + t4 * t7, -t4, 0, 0;];
	Ja_transl = t8;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:38:50
	% EndTime: 2019-12-31 17:38:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (20->9), mult. (43->12), div. (0->0), fcn. (54->6), ass. (0->10)
	t1 = sin(pkin(8));
	t3 = cos(pkin(8));
	t9 = t1 * r_i_i_C(1) + t3 * r_i_i_C(2) + qJ(3);
	t8 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2) + pkin(2) + pkin(3);
	t5 = sin(qJ(2));
	t6 = cos(qJ(2));
	t7 = -t8 * t5 + t9 * t6;
	t4 = cos(pkin(7));
	t2 = sin(pkin(7));
	t10 = [0, t7 * t4, t4 * t5, -t2, 0; 0, t7 * t2, t2 * t5, t4, 0; 1, t9 * t5 + t8 * t6, -t6, 0, 0;];
	Ja_transl = t10;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:38:50
	% EndTime: 2019-12-31 17:38:50
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (51->23), mult. (124->36), div. (0->0), fcn. (154->8), ass. (0->17)
	t16 = sin(qJ(2));
	t18 = cos(qJ(2));
	t15 = sin(qJ(5));
	t17 = cos(qJ(5));
	t21 = r_i_i_C(1) * t17 - r_i_i_C(2) * t15 + pkin(4);
	t11 = sin(pkin(8));
	t13 = cos(pkin(8));
	t22 = t11 * t18 - t16 * t13;
	t26 = pkin(2) + pkin(3);
	t28 = qJ(3) * t18 - t26 * t16 + t22 * t21;
	t27 = t16 * t11 + t13 * t18;
	t25 = pkin(6) + r_i_i_C(3);
	t14 = cos(pkin(7));
	t12 = sin(pkin(7));
	t3 = t27 * t14;
	t1 = t27 * t12;
	t2 = [0, t28 * t14 - t25 * t3, t14 * t16, -t12, (-t12 * t17 - t15 * t3) * r_i_i_C(1) + (t12 * t15 - t17 * t3) * r_i_i_C(2); 0, -t25 * t1 + t28 * t12, t12 * t16, t14, (-t1 * t15 + t14 * t17) * r_i_i_C(1) + (-t1 * t17 - t14 * t15) * r_i_i_C(2); 1, t16 * qJ(3) + t26 * t18 + t21 * t27 + t22 * t25, -t18, 0, -(-r_i_i_C(1) * t15 - r_i_i_C(2) * t17) * t22;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,5);
end