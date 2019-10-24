% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:29
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRPP3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPP3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
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
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (19->12), mult. (51->25), div. (0->0), fcn. (55->6), ass. (0->12)
	t3 = sin(qJ(3));
	t6 = cos(qJ(2));
	t11 = t3 * t6;
	t5 = cos(qJ(3));
	t10 = t5 * t6;
	t9 = pkin(6) + r_i_i_C(3);
	t8 = t5 * r_i_i_C(1) - t3 * r_i_i_C(2) + pkin(2);
	t4 = sin(qJ(2));
	t7 = -t8 * t4 + t9 * t6;
	t2 = cos(pkin(7));
	t1 = sin(pkin(7));
	t12 = [0, t7 * t2, (t1 * t5 - t2 * t11) * r_i_i_C(1) + (-t1 * t3 - t2 * t10) * r_i_i_C(2), 0, 0; 0, t7 * t1, (-t1 * t11 - t2 * t5) * r_i_i_C(1) + (-t1 * t10 + t2 * t3) * r_i_i_C(2), 0, 0; 1, t9 * t4 + t8 * t6, (-r_i_i_C(1) * t3 - r_i_i_C(2) * t5) * t4, 0, 0;];
	Ja_transl = t12;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (48->16), mult. (129->30), div. (0->0), fcn. (150->8), ass. (0->18)
	t10 = sin(qJ(3));
	t12 = cos(qJ(3));
	t6 = sin(pkin(8));
	t8 = cos(pkin(8));
	t16 = r_i_i_C(1) * t8 - r_i_i_C(2) * t6 + pkin(3);
	t17 = r_i_i_C(3) + qJ(4);
	t20 = t17 * t10 + t16 * t12 + pkin(2);
	t13 = cos(qJ(2));
	t19 = t10 * t13;
	t18 = t12 * t13;
	t15 = t6 * r_i_i_C(1) + t8 * r_i_i_C(2) + pkin(6);
	t11 = sin(qJ(2));
	t14 = -t20 * t11 + t15 * t13;
	t9 = cos(pkin(7));
	t7 = sin(pkin(7));
	t3 = -t7 * t12 + t9 * t19;
	t1 = t9 * t12 + t7 * t19;
	t2 = [0, t14 * t9, t17 * (t7 * t10 + t9 * t18) - t16 * t3, t3, 0; 0, t14 * t7, t17 * (-t9 * t10 + t7 * t18) - t16 * t1, t1, 0; 1, t15 * t11 + t20 * t13, (-t16 * t10 + t17 * t12) * t11, t11 * t10, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:29:31
	% EndTime: 2019-10-24 10:29:31
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (75->26), mult. (202->45), div. (0->0), fcn. (241->8), ass. (0->24)
	t12 = sin(pkin(8));
	t14 = cos(pkin(8));
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t18 = cos(qJ(3));
	t29 = t17 * t18;
	t24 = t12 * t29 + t14 * t19;
	t26 = r_i_i_C(3) + qJ(5);
	t32 = pkin(4) + r_i_i_C(1);
	t16 = sin(qJ(3));
	t27 = r_i_i_C(2) + qJ(4);
	t33 = pkin(3) * t18 + t27 * t16 + pkin(2);
	t34 = pkin(6) * t19 - t33 * t17 - t26 * t24 + t32 * (t12 * t19 - t14 * t29);
	t31 = t16 * t19;
	t30 = t17 * t14;
	t28 = t18 * t19;
	t23 = -t26 * t12 - t32 * t14 - pkin(3);
	t15 = cos(pkin(7));
	t13 = sin(pkin(7));
	t8 = t13 * t16 + t15 * t28;
	t7 = -t13 * t18 + t15 * t31;
	t6 = t13 * t28 - t15 * t16;
	t5 = t13 * t31 + t15 * t18;
	t1 = [0, t34 * t15, t23 * t7 + t27 * t8, t7, t12 * t8 - t15 * t30; 0, t34 * t13, t23 * t5 + t27 * t6, t5, t12 * t6 - t13 * t30; 1, t17 * pkin(6) + t26 * (t12 * t28 - t30) + t32 * (t17 * t12 + t14 * t28) + t33 * t19, (t23 * t16 + t27 * t18) * t17, t17 * t16, t24;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end