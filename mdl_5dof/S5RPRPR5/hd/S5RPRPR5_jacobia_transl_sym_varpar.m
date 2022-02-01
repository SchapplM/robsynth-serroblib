% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRPR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
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
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (24->17), mult. (48->25), div. (0->0), fcn. (58->6), ass. (0->16)
	t8 = sin(qJ(3));
	t9 = sin(qJ(1));
	t16 = t9 * t8;
	t11 = cos(qJ(1));
	t15 = t11 * t8;
	t10 = cos(qJ(3));
	t14 = t9 * t10;
	t13 = t11 * t10;
	t6 = sin(pkin(8));
	t7 = cos(pkin(8));
	t12 = pkin(2) * t7 + pkin(1) + (pkin(6) + r_i_i_C(3)) * t6;
	t4 = t7 * t13 + t16;
	t3 = -t7 * t15 + t14;
	t2 = -t7 * t14 + t15;
	t1 = t7 * t16 + t13;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t11 * qJ(2) - t12 * t9, t9, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t9 * qJ(2) + t12 * t11, -t11, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0; 0, 0, (-r_i_i_C(1) * t8 - r_i_i_C(2) * t10) * t6, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (54->27), mult. (68->37), div. (0->0), fcn. (81->8), ass. (0->19)
	t12 = sin(qJ(3));
	t19 = t12 * pkin(3);
	t11 = cos(pkin(8));
	t13 = sin(qJ(1));
	t18 = t11 * t13;
	t15 = cos(qJ(1));
	t17 = t11 * t15;
	t10 = sin(pkin(8));
	t14 = cos(qJ(3));
	t16 = pkin(1) + (t14 * pkin(3) + pkin(2)) * t11 + (r_i_i_C(3) + qJ(4) + pkin(6)) * t10;
	t9 = qJ(3) + pkin(9);
	t8 = cos(t9);
	t7 = sin(t9);
	t6 = qJ(2) + t19;
	t4 = t13 * t7 + t8 * t17;
	t3 = t13 * t8 - t7 * t17;
	t2 = t15 * t7 - t8 * t18;
	t1 = t15 * t8 + t7 * t18;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t16 * t13 + t6 * t15, t13, t3 * r_i_i_C(1) - t4 * r_i_i_C(2) + (-t12 * t17 + t13 * t14) * pkin(3), t10 * t15, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t6 * t13 + t16 * t15, -t15, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + (-t12 * t18 - t14 * t15) * pkin(3), t10 * t13, 0; 0, 0, (-r_i_i_C(1) * t7 - r_i_i_C(2) * t8 - t19) * t10, -t11, 0;];
	Ja_transl = t5;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (119->31), mult. (101->38), div. (0->0), fcn. (118->10), ass. (0->23)
	t18 = qJ(3) + pkin(9);
	t11 = pkin(4) * cos(t18) + cos(qJ(3)) * pkin(3);
	t19 = sin(pkin(8));
	t20 = cos(pkin(8));
	t32 = (r_i_i_C(3) + pkin(7) + qJ(4) + pkin(6)) * t19 + t20 * (pkin(2) + t11) + pkin(1);
	t15 = qJ(5) + t18;
	t14 = cos(t15);
	t22 = cos(qJ(1));
	t13 = sin(t15);
	t21 = sin(qJ(1));
	t27 = t21 * t13;
	t5 = t14 * t22 + t20 * t27;
	t26 = t21 * t14;
	t6 = t13 * t22 - t20 * t26;
	t31 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t28 = t20 * t22;
	t7 = -t13 * t28 + t26;
	t8 = t14 * t28 + t27;
	t30 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t10 = pkin(4) * sin(t18) + sin(qJ(3)) * pkin(3);
	t25 = qJ(2) + t10;
	t23 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t32 * t21 + t25 * t22, t21, -t10 * t28 + t21 * t11 + t30, t22 * t19, t30; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t25 * t21 + t32 * t22, -t22, -t21 * t20 * t10 - t11 * t22 + t31, t21 * t19, t31; 0, 0, (-t10 + t23) * t19, -t20, t23 * t19;];
	Ja_transl = t1;
end