% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPR8
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
% Datum: 2019-12-29 15:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRPR8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPR8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:16
	% EndTime: 2019-12-29 15:40:16
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:21
	% EndTime: 2019-12-29 15:40:21
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:16
	% EndTime: 2019-12-29 15:40:16
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (3->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t3 = sin(qJ(2));
	t4 = cos(qJ(2));
	t5 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t1 = [0, t5 * cos(pkin(8)), 0, 0, 0; 0, t5 * sin(pkin(8)), 0, 0, 0; 1, r_i_i_C(1) * t4 - r_i_i_C(2) * t3, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:16
	% EndTime: 2019-12-29 15:40:16
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (21->5), mult. (25->10), div. (0->0), fcn. (25->6), ass. (0->9)
	t4 = qJ(2) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = -sin(qJ(2)) * pkin(2) + t9;
	t6 = cos(pkin(8));
	t5 = sin(pkin(8));
	t1 = [0, t8 * t6, t9 * t6, 0, 0; 0, t8 * t5, t9 * t5, 0, 0; 1, cos(qJ(2)) * pkin(2) + t10, t10, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:21
	% EndTime: 2019-12-29 15:40:21
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (46->9), mult. (33->12), div. (0->0), fcn. (35->8), ass. (0->10)
	t8 = qJ(2) + qJ(3);
	t6 = pkin(9) + t8;
	t3 = sin(t6);
	t4 = cos(t6);
	t11 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4 - pkin(3) * sin(t8);
	t14 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + pkin(3) * cos(t8);
	t12 = -sin(qJ(2)) * pkin(2) + t11;
	t10 = cos(pkin(8));
	t9 = sin(pkin(8));
	t1 = [0, t12 * t10, t11 * t10, t9, 0; 0, t12 * t9, t11 * t9, -t10, 0; 1, cos(qJ(2)) * pkin(2) + t14, t14, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:40:21
	% EndTime: 2019-12-29 15:40:21
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (119->23), mult. (95->32), div. (0->0), fcn. (101->10), ass. (0->22)
	t17 = qJ(2) + qJ(3);
	t15 = pkin(9) + t17;
	t12 = sin(t15);
	t13 = cos(t15);
	t20 = sin(qJ(5));
	t35 = r_i_i_C(2) * t20;
	t40 = pkin(7) + r_i_i_C(3);
	t41 = t12 * t35 + t13 * t40;
	t21 = cos(qJ(5));
	t39 = -r_i_i_C(1) * t21 - pkin(4);
	t23 = t39 * t12 - pkin(3) * sin(t17);
	t18 = sin(pkin(8));
	t32 = t18 * t20;
	t31 = t18 * t21;
	t19 = cos(pkin(8));
	t30 = t19 * t20;
	t29 = t19 * t21;
	t28 = t41 * t18;
	t27 = t41 * t19;
	t24 = -sin(qJ(2)) * pkin(2) + t23;
	t22 = pkin(3) * cos(t17) + t40 * t12 + (-t35 - t39) * t13;
	t1 = [0, t24 * t19 + t27, t23 * t19 + t27, t18, (-t13 * t30 + t31) * r_i_i_C(1) + (-t13 * t29 - t32) * r_i_i_C(2); 0, t24 * t18 + t28, t23 * t18 + t28, -t19, (-t13 * t32 - t29) * r_i_i_C(1) + (-t13 * t31 + t30) * r_i_i_C(2); 1, cos(qJ(2)) * pkin(2) + t22, t22, 0, (-r_i_i_C(1) * t20 - r_i_i_C(2) * t21) * t12;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end