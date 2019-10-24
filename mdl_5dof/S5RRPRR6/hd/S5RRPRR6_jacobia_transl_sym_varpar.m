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
% Datum: 2019-10-24 10:49
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:49:50
	% EndTime: 2019-10-24 10:49:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:50
	% EndTime: 2019-10-24 10:49:50
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
	% StartTime: 2019-10-24 10:49:50
	% EndTime: 2019-10-24 10:49:50
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
	% StartTime: 2019-10-24 10:49:50
	% EndTime: 2019-10-24 10:49:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (40->8), mult. (30->8), div. (0->0), fcn. (32->6), ass. (0->8)
	t15 = r_i_i_C(3) + qJ(3);
	t14 = r_i_i_C(2) * sin(pkin(9)) - r_i_i_C(1) * cos(pkin(9)) - pkin(2);
	t7 = qJ(1) + qJ(2);
	t5 = sin(t7);
	t6 = cos(t7);
	t11 = t14 * t5 + t15 * t6;
	t10 = t14 * t6 - t15 * t5;
	t1 = [0, 0, 0, 0, 0; -cos(qJ(1)) * pkin(1) + t10, t10, t6, 0, 0; -sin(qJ(1)) * pkin(1) + t11, t11, t5, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:50
	% EndTime: 2019-10-24 10:49:50
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (87->19), mult. (86->29), div. (0->0), fcn. (100->8), ass. (0->17)
	t14 = cos(pkin(9));
	t15 = sin(qJ(4));
	t21 = t14 * t15;
	t16 = cos(qJ(4));
	t20 = t14 * t16;
	t13 = sin(pkin(9));
	t19 = -pkin(3) * t14 - pkin(2) + (-pkin(7) - r_i_i_C(3)) * t13;
	t12 = qJ(1) + qJ(2);
	t10 = sin(t12);
	t11 = cos(t12);
	t5 = t10 * t21 + t11 * t16;
	t6 = t10 * t20 - t11 * t15;
	t18 = -t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t11 * qJ(3) + t19 * t10;
	t7 = -t10 * t16 + t11 * t21;
	t8 = -t10 * t15 - t11 * t20;
	t17 = t8 * r_i_i_C(1) + t7 * r_i_i_C(2) - t10 * qJ(3) + t19 * t11;
	t1 = [0, 0, 0, (-r_i_i_C(1) * t15 - r_i_i_C(2) * t16) * t13, 0; -cos(qJ(1)) * pkin(1) + t17, t17, t11, r_i_i_C(1) * t5 + r_i_i_C(2) * t6, 0; -sin(qJ(1)) * pkin(1) + t18, t18, t10, -r_i_i_C(1) * t7 + r_i_i_C(2) * t8, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:49:50
	% EndTime: 2019-10-24 10:49:50
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (163->28), mult. (130->39), div. (0->0), fcn. (148->10), ass. (0->26)
	t24 = sin(pkin(9));
	t25 = cos(pkin(9));
	t27 = cos(qJ(4));
	t41 = -(t27 * pkin(4) + pkin(3)) * t25 - pkin(2) + (-pkin(8) - pkin(7) - r_i_i_C(3)) * t24;
	t26 = sin(qJ(4));
	t37 = pkin(4) * t26;
	t40 = qJ(3) + t37;
	t22 = qJ(4) + qJ(5);
	t18 = sin(t22);
	t20 = cos(t22);
	t23 = qJ(1) + qJ(2);
	t21 = cos(t23);
	t19 = sin(t23);
	t36 = t19 * t25;
	t10 = -t21 * t18 + t20 * t36;
	t9 = t18 * t36 + t21 * t20;
	t39 = t9 * r_i_i_C(1) + t10 * r_i_i_C(2);
	t35 = t21 * t25;
	t11 = t18 * t35 - t19 * t20;
	t12 = -t19 * t18 - t20 * t35;
	t38 = -t11 * r_i_i_C(1) + t12 * r_i_i_C(2);
	t33 = t25 * t26;
	t32 = -r_i_i_C(1) * t18 - r_i_i_C(2) * t20;
	t30 = -t10 * r_i_i_C(1) + t9 * r_i_i_C(2) + t41 * t19 + t40 * t21;
	t29 = t12 * r_i_i_C(1) + t11 * r_i_i_C(2) - t40 * t19 + t41 * t21;
	t1 = [0, 0, 0, (t32 - t37) * t24, t32 * t24; -cos(qJ(1)) * pkin(1) + t29, t29, t21, (t19 * t33 + t21 * t27) * pkin(4) + t39, t39; -sin(qJ(1)) * pkin(1) + t30, t30, t19, (t19 * t27 - t21 * t33) * pkin(4) + t38, t38;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end