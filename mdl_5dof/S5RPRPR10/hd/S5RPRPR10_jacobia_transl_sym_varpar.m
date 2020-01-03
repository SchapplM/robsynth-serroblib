% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRPR10_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR10_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR10_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->5), mult. (8->4), div. (0->0), fcn. (10->2), ass. (0->5)
	t4 = pkin(1) + r_i_i_C(1);
	t3 = r_i_i_C(3) + qJ(2);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t4 * t1 + t3 * t2, t1, 0, 0, 0; t3 * t1 + t4 * t2, -t2, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (19->10), mult. (30->12), div. (0->0), fcn. (40->4), ass. (0->10)
	t12 = pkin(1) + pkin(2);
	t11 = cos(qJ(3));
	t10 = sin(qJ(3));
	t6 = sin(qJ(1));
	t7 = cos(qJ(1));
	t1 = -t6 * t10 - t7 * t11;
	t2 = t7 * t10 - t6 * t11;
	t9 = -t2 * r_i_i_C(1) + t1 * r_i_i_C(2);
	t8 = t1 * r_i_i_C(1) + t2 * r_i_i_C(2);
	t3 = [t7 * qJ(2) - t12 * t6 - t9, t6, t9, 0, 0; t6 * qJ(2) + t12 * t7 - t8, -t7, t8, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (43->17), mult. (44->20), div. (0->0), fcn. (54->6), ass. (0->14)
	t9 = cos(qJ(3));
	t17 = t9 * pkin(3) + pkin(1) + pkin(2);
	t16 = qJ(3) + pkin(8);
	t7 = sin(qJ(3));
	t15 = pkin(3) * t7 + qJ(2);
	t14 = cos(t16);
	t13 = sin(t16);
	t10 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = -t10 * t14 - t8 * t13;
	t2 = t10 * t13 - t8 * t14;
	t12 = -t2 * r_i_i_C(1) + t1 * r_i_i_C(2);
	t11 = t1 * r_i_i_C(1) + t2 * r_i_i_C(2);
	t3 = [t15 * t10 - t17 * t8 - t12, t8, (-t10 * t7 + t8 * t9) * pkin(3) + t12, 0, 0; t17 * t10 + t15 * t8 - t11, -t10, (-t10 * t9 - t7 * t8) * pkin(3) + t11, 0, 0; 0, 0, 0, -1, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (106->22), mult. (106->27), div. (0->0), fcn. (132->8), ass. (0->20)
	t10 = cos(qJ(3));
	t23 = pkin(3) * t10 + pkin(1) + pkin(2);
	t11 = cos(qJ(1));
	t18 = qJ(3) + pkin(8);
	t15 = sin(t18);
	t16 = cos(t18);
	t19 = sin(qJ(1));
	t1 = -t11 * t16 - t19 * t15;
	t7 = sin(qJ(5));
	t9 = cos(qJ(5));
	t14 = -r_i_i_C(1) * t9 + r_i_i_C(2) * t7;
	t12 = pkin(4) - t14;
	t2 = t11 * t15 - t19 * t16;
	t20 = pkin(7) + r_i_i_C(3);
	t22 = -t20 * t1 - t12 * t2;
	t21 = t12 * t1 - t20 * t2;
	t8 = sin(qJ(3));
	t17 = t19 * t8;
	t13 = r_i_i_C(1) * t7 + r_i_i_C(2) * t9;
	t3 = [(pkin(3) * t8 + qJ(2)) * t11 - t23 * t19 - t22, t19, (t19 * t10 - t11 * t8) * pkin(3) + t22, 0, t13 * t1; pkin(3) * t17 + t19 * qJ(2) + t23 * t11 - t21, -t11, (-t10 * t11 - t17) * pkin(3) + t21, 0, t13 * t2; 0, 0, 0, -1, t14;];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,5);
end