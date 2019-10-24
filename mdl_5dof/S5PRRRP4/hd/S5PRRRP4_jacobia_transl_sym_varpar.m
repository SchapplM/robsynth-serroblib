% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:33
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRRP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:33:46
	% EndTime: 2019-10-24 10:33:46
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:33:46
	% EndTime: 2019-10-24 10:33:46
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:33:46
	% EndTime: 2019-10-24 10:33:46
	% DurationCPUTime: 0.06s
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
	% StartTime: 2019-10-24 10:33:46
	% EndTime: 2019-10-24 10:33:46
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (21->5), mult. (25->10), div. (0->0), fcn. (25->6), ass. (0->9)
	t4 = qJ(2) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = r_i_i_C(1) * t3 - r_i_i_C(2) * t2;
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
	% StartTime: 2019-10-24 10:33:46
	% EndTime: 2019-10-24 10:33:46
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
	t14 = sin(pkin(8));
	t28 = t14 * t16;
	t27 = t14 * t18;
	t15 = cos(pkin(8));
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
	% StartTime: 2019-10-24 10:33:46
	% EndTime: 2019-10-24 10:33:46
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (115->20), mult. (146->30), div. (0->0), fcn. (159->8), ass. (0->24)
	t22 = sin(qJ(4));
	t24 = cos(qJ(4));
	t28 = r_i_i_C(3) + qJ(5);
	t43 = pkin(4) + r_i_i_C(1);
	t45 = t28 * t22 + t24 * t43 + pkin(3);
	t42 = pkin(7) + r_i_i_C(2);
	t19 = qJ(2) + qJ(3);
	t18 = cos(t19);
	t40 = t18 * t42;
	t20 = sin(pkin(8));
	t38 = t20 * t40;
	t21 = cos(pkin(8));
	t37 = t21 * t40;
	t32 = t20 * t22;
	t31 = t20 * t24;
	t30 = t21 * t22;
	t29 = t21 * t24;
	t17 = sin(t19);
	t27 = t42 * t17 + t45 * t18;
	t26 = t45 * t17;
	t25 = -sin(qJ(2)) * pkin(2) - t26;
	t3 = t18 * t30 - t31;
	t1 = t18 * t32 + t29;
	t2 = [0, t25 * t21 + t37, -t21 * t26 + t37, t28 * (t18 * t29 + t32) - t43 * t3, t3; 0, t25 * t20 + t38, -t20 * t26 + t38, t28 * (t18 * t31 - t30) - t43 * t1, t1; 1, cos(qJ(2)) * pkin(2) + t27, t27, (-t22 * t43 + t28 * t24) * t17, t17 * t22;];
	Ja_transl = t2;
else
	Ja_transl=NaN(3,5);
end