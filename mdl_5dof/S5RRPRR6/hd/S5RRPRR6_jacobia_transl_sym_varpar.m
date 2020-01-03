% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPRR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t6 = qJ(1) + qJ(2);
	t4 = sin(t6);
	t5 = cos(t6);
	t8 = t4 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t7 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t1 = [0, 0, 0, 0, 0; cos(qJ(1)) * pkin(1) + t7, t7, 0, 0, 0; sin(qJ(1)) * pkin(1) + t8, t8, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (42->10), mult. (30->8), div. (0->0), fcn. (32->6), ass. (0->8)
	t17 = r_i_i_C(3) + qJ(3);
	t16 = -r_i_i_C(2) * sin(pkin(9)) + r_i_i_C(1) * cos(pkin(9)) + pkin(2);
	t9 = qJ(1) + qJ(2);
	t7 = sin(t9);
	t8 = cos(t9);
	t13 = t16 * t8 + t17 * t7;
	t12 = t16 * t7 - t17 * t8;
	t1 = [0, 0, 0, 0, 0; cos(qJ(1)) * pkin(1) + t13, t13, -t8, 0, 0; sin(qJ(1)) * pkin(1) + t12, t12, -t7, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (89->21), mult. (86->29), div. (0->0), fcn. (100->8), ass. (0->17)
	t21 = sin(pkin(9));
	t22 = cos(pkin(9));
	t34 = t21 * (pkin(7) + r_i_i_C(3)) + pkin(3) * t22 + pkin(2);
	t23 = sin(qJ(4));
	t28 = t22 * t23;
	t24 = cos(qJ(4));
	t27 = t22 * t24;
	t20 = qJ(1) + qJ(2);
	t18 = sin(t20);
	t19 = cos(t20);
	t7 = -t18 * t24 + t19 * t28;
	t8 = t18 * t23 + t19 * t27;
	t26 = t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t18 * qJ(3) + t34 * t19;
	t5 = -t18 * t28 - t19 * t24;
	t6 = t18 * t27 - t19 * t23;
	t25 = t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t19 * qJ(3) + t34 * t18;
	t1 = [0, 0, 0, (-r_i_i_C(1) * t23 - r_i_i_C(2) * t24) * t21, 0; cos(qJ(1)) * pkin(1) + t26, t26, -t19, t5 * r_i_i_C(1) - t6 * r_i_i_C(2), 0; sin(qJ(1)) * pkin(1) + t25, t25, -t18, t7 * r_i_i_C(1) + t8 * r_i_i_C(2), 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:06:32
	% EndTime: 2020-01-03 12:06:32
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (165->31), mult. (130->40), div. (0->0), fcn. (148->10), ass. (0->27)
	t30 = sin(qJ(4));
	t41 = pkin(4) * t30;
	t45 = qJ(3) + t41;
	t28 = sin(pkin(9));
	t44 = pkin(2) + (pkin(8) + pkin(7) + r_i_i_C(3)) * t28;
	t26 = qJ(4) + qJ(5);
	t22 = sin(t26);
	t24 = cos(t26);
	t27 = qJ(1) + qJ(2);
	t25 = cos(t27);
	t23 = sin(t27);
	t29 = cos(pkin(9));
	t39 = t23 * t29;
	t10 = -t22 * t25 + t24 * t39;
	t9 = -t22 * t39 - t24 * t25;
	t43 = t9 * r_i_i_C(1) - t10 * r_i_i_C(2);
	t38 = t25 * t29;
	t11 = t22 * t38 - t23 * t24;
	t12 = t22 * t23 + t24 * t38;
	t42 = t11 * r_i_i_C(1) + t12 * r_i_i_C(2);
	t36 = t29 * t30;
	t35 = -r_i_i_C(1) * t22 - r_i_i_C(2) * t24;
	t31 = cos(qJ(4));
	t21 = pkin(4) * t31 + pkin(3);
	t34 = t12 * r_i_i_C(1) - t11 * r_i_i_C(2) + t21 * t38 + t45 * t23 + t44 * t25;
	t33 = t10 * r_i_i_C(1) + t9 * r_i_i_C(2) + t21 * t39 + t44 * t23 - t45 * t25;
	t1 = [0, 0, 0, (t35 - t41) * t28, t35 * t28; cos(qJ(1)) * pkin(1) + t34, t34, -t25, (-t23 * t36 - t25 * t31) * pkin(4) + t43, t43; sin(qJ(1)) * pkin(1) + t33, t33, -t23, (-t23 * t31 + t25 * t36) * pkin(4) + t42, t42;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end