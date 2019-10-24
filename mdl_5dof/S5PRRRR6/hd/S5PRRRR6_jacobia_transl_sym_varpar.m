% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:37
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRRR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:19
	% EndTime: 2019-10-24 10:37:19
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:19
	% EndTime: 2019-10-24 10:37:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:19
	% EndTime: 2019-10-24 10:37:19
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t3 = sin(qJ(2));
	t4 = cos(qJ(2));
	t5 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t1 = [0, t5 * cos(pkin(9)), 0, 0, 0; 0, t5 * sin(pkin(9)), 0, 0, 0; 1, r_i_i_C(1) * t4 - r_i_i_C(2) * t3, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:19
	% EndTime: 2019-10-24 10:37:19
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (21->5), mult. (25->10), div. (0->0), fcn. (25->6), ass. (0->9)
	t4 = qJ(2) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = r_i_i_C(1) * t3 - r_i_i_C(2) * t2;
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = -sin(qJ(2)) * pkin(2) + t9;
	t6 = cos(pkin(9));
	t5 = sin(pkin(9));
	t1 = [0, t8 * t6, t9 * t6, 0, 0; 0, t8 * t5, t9 * t5, 0, 0; 1, cos(qJ(2)) * pkin(2) + t10, t10, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:19
	% EndTime: 2019-10-24 10:37:19
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (70->19), mult. (87->30), div. (0->0), fcn. (91->8), ass. (0->21)
	t13 = qJ(2) + qJ(3);
	t11 = sin(t13);
	t12 = cos(t13);
	t16 = sin(qJ(4));
	t31 = r_i_i_C(2) * t16;
	t35 = pkin(7) + r_i_i_C(3);
	t36 = t11 * t31 + t12 * t35;
	t18 = cos(qJ(4));
	t34 = -r_i_i_C(1) * t18 - pkin(3);
	t14 = sin(pkin(9));
	t28 = t14 * t16;
	t27 = t14 * t18;
	t15 = cos(pkin(9));
	t26 = t15 * t16;
	t25 = t15 * t18;
	t24 = t36 * t14;
	t23 = t36 * t15;
	t21 = t34 * t11;
	t20 = t35 * t11 + (-t31 - t34) * t12;
	t19 = -sin(qJ(2)) * pkin(2) + t21;
	t1 = [0, t19 * t15 + t23, t15 * t21 + t23, (-t12 * t26 + t27) * r_i_i_C(1) + (-t12 * t25 - t28) * r_i_i_C(2), 0; 0, t19 * t14 + t24, t14 * t21 + t24, (-t12 * t28 - t25) * r_i_i_C(1) + (-t12 * t27 + t26) * r_i_i_C(2), 0; 1, cos(qJ(2)) * pkin(2) + t20, t20, (-r_i_i_C(1) * t16 - r_i_i_C(2) * t18) * t11, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:37:19
	% EndTime: 2019-10-24 10:37:19
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (135->29), mult. (125->44), div. (0->0), fcn. (133->10), ass. (0->26)
	t17 = qJ(4) + qJ(5);
	t15 = cos(t17);
	t23 = cos(qJ(4));
	t39 = -pkin(4) * t23 - r_i_i_C(1) * t15 - pkin(3);
	t13 = sin(t17);
	t20 = cos(pkin(9));
	t18 = qJ(2) + qJ(3);
	t16 = cos(t18);
	t19 = sin(pkin(9));
	t32 = t16 * t19;
	t38 = (-t13 * t32 - t15 * t20) * r_i_i_C(1) + (t13 * t20 - t15 * t32) * r_i_i_C(2);
	t31 = t16 * t20;
	t37 = (-t13 * t31 + t15 * t19) * r_i_i_C(1) + (-t13 * t19 - t15 * t31) * r_i_i_C(2);
	t14 = sin(t18);
	t34 = r_i_i_C(2) * t13;
	t29 = t14 * t34;
	t36 = r_i_i_C(3) * t32 + t19 * t29;
	t33 = r_i_i_C(3) * t31 + t20 * t29;
	t21 = sin(qJ(4));
	t30 = t16 * t21;
	t28 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t15;
	t24 = -pkin(8) - pkin(7);
	t27 = t39 * t14 - t16 * t24;
	t26 = (r_i_i_C(3) - t24) * t14 + (-t34 - t39) * t16;
	t25 = -sin(qJ(2)) * pkin(2) + t27;
	t1 = [0, t25 * t20 + t33, t27 * t20 + t33, (t19 * t23 - t20 * t30) * pkin(4) + t37, t37; 0, t25 * t19 + t36, t27 * t19 + t36, (-t19 * t30 - t20 * t23) * pkin(4) + t38, t38; 1, cos(qJ(2)) * pkin(2) + t26, t26, (-pkin(4) * t21 + t28) * t14, t28 * t14;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end