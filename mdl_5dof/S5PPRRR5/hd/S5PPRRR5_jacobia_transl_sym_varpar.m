% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRRR5
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PPRRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:35:59
	% EndTime: 2019-12-31 17:35:59
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:35:59
	% EndTime: 2019-12-31 17:35:59
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:35:59
	% EndTime: 2019-12-31 17:35:59
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (2->2), ass. (0->1)
	t1 = [0, sin(pkin(8)), 0, 0, 0; 0, -cos(pkin(8)), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:35:59
	% EndTime: 2019-12-31 17:35:59
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (7->5), mult. (12->8), div. (0->0), fcn. (18->4), ass. (0->7)
	t6 = cos(qJ(3));
	t5 = sin(qJ(3));
	t4 = cos(pkin(8));
	t3 = sin(pkin(8));
	t2 = -t3 * t6 + t4 * t5;
	t1 = -t3 * t5 - t4 * t6;
	t7 = [0, t3, -t2 * r_i_i_C(1) + t1 * r_i_i_C(2), 0, 0; 0, -t4, t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:35:59
	% EndTime: 2019-12-31 17:35:59
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (33->10), mult. (32->14), div. (0->0), fcn. (42->6), ass. (0->12)
	t10 = sin(pkin(8));
	t11 = cos(pkin(8));
	t9 = qJ(3) + qJ(4);
	t7 = sin(t9);
	t8 = cos(t9);
	t5 = -t10 * t7 - t11 * t8;
	t6 = -t10 * t8 + t11 * t7;
	t15 = t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t14 = -t6 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t13 = cos(qJ(3));
	t12 = sin(qJ(3));
	t1 = [0, t10, (t10 * t13 - t11 * t12) * pkin(3) + t14, t14, 0; 0, -t11, (-t10 * t12 - t11 * t13) * pkin(3) + t15, t15, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:35:59
	% EndTime: 2019-12-31 17:35:59
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (96->14), mult. (94->20), div. (0->0), fcn. (120->8), ass. (0->18)
	t28 = -pkin(7) - r_i_i_C(3);
	t13 = sin(qJ(5));
	t15 = cos(qJ(5));
	t27 = -t15 * r_i_i_C(1) + t13 * r_i_i_C(2);
	t26 = pkin(4) - t27;
	t23 = sin(pkin(8));
	t22 = qJ(3) + qJ(4);
	t21 = cos(t22);
	t20 = sin(t22);
	t19 = r_i_i_C(1) * t13 + r_i_i_C(2) * t15;
	t12 = cos(pkin(8));
	t7 = -t12 * t21 - t23 * t20;
	t8 = t12 * t20 - t23 * t21;
	t18 = t26 * t7 + t28 * t8;
	t17 = -t26 * t8 + t28 * t7;
	t16 = cos(qJ(3));
	t14 = sin(qJ(3));
	t1 = [0, t23, (-t12 * t14 + t23 * t16) * pkin(3) + t17, t17, t19 * t7; 0, -t12, (-t12 * t16 - t23 * t14) * pkin(3) + t18, t18, t19 * t8; 1, 0, 0, 0, t27;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end