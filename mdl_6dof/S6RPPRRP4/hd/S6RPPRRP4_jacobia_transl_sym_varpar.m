% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPPRRP4
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:52
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->5), mult. (8->4), div. (0->0), fcn. (10->2), ass. (0->5)
	t4 = pkin(1) + r_i_i_C(1);
	t3 = r_i_i_C(3) + qJ(2);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t1 * t4 + t2 * t3, t1, 0, 0, 0, 0; t1 * t3 + t2 * t4, -t2, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (13->10), mult. (18->12), div. (0->0), fcn. (24->4), ass. (0->8)
	t7 = pkin(1) + pkin(2);
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t4 = cos(pkin(9));
	t3 = sin(pkin(9));
	t2 = t3 * t5 + t4 * t6;
	t1 = t3 * t6 - t4 * t5;
	t8 = [t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + t6 * qJ(2) - t5 * t7, t5, 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + t5 * qJ(2) + t6 * t7, -t6, 0, 0, 0, 0; 0, 0, -1, 0, 0, 0;];
	Ja_transl = t8;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (32->14), mult. (58->18), div. (0->0), fcn. (74->6), ass. (0->14)
	t16 = pkin(1) + pkin(2);
	t15 = pkin(7) + r_i_i_C(3);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t12 = cos(pkin(9));
	t11 = sin(pkin(9));
	t6 = sin(qJ(4));
	t7 = cos(qJ(4));
	t10 = -t7 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t9 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7;
	t8 = pkin(3) - t10;
	t2 = t14 * t11 - t13 * t12;
	t1 = -t13 * t11 - t14 * t12;
	t3 = [t14 * qJ(2) + t15 * t1 - t16 * t13 + t8 * t2, t13, 0, t9 * t1, 0, 0; t13 * qJ(2) - t8 * t1 + t16 * t14 + t15 * t2, -t14, 0, t9 * t2, 0, 0; 0, 0, -1, t10, 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (78->28), mult. (161->40), div. (0->0), fcn. (207->8), ass. (0->22)
	t11 = cos(qJ(4));
	t10 = cos(qJ(5));
	t8 = sin(qJ(5));
	t13 = r_i_i_C(1) * t10 - r_i_i_C(2) * t8 + pkin(4);
	t22 = pkin(8) + r_i_i_C(3);
	t9 = sin(qJ(4));
	t15 = t22 * t9;
	t25 = -t13 * t11 - t15;
	t24 = pkin(1) + pkin(2);
	t21 = t11 * t8;
	t20 = cos(qJ(1));
	t19 = sin(qJ(1));
	t18 = t10 * t11;
	t17 = cos(pkin(9));
	t16 = sin(pkin(9));
	t14 = t8 * r_i_i_C(1) + t10 * r_i_i_C(2);
	t12 = -t22 * t11 + t13 * t9;
	t4 = t20 * t16 - t19 * t17;
	t3 = -t19 * t16 - t20 * t17;
	t2 = -t3 * t18 + t4 * t8;
	t1 = t4 * t10 + t3 * t21;
	t5 = [t20 * qJ(2) + (pkin(7) + t14) * t3 + (pkin(3) - t25) * t4 - t24 * t19, t19, 0, t12 * t3, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t19 * qJ(2) + t4 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + (-pkin(4) * t11 - pkin(3) - t15) * t3 + t24 * t20, -t20, 0, t12 * t4, (-t3 * t10 + t4 * t21) * r_i_i_C(1) + (t4 * t18 + t3 * t8) * r_i_i_C(2), 0; 0, 0, -1, t25, t14 * t9, 0;];
	Ja_transl = t5;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:52:56
	% EndTime: 2019-10-09 23:52:56
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (129->32), mult. (268->43), div. (0->0), fcn. (351->8), ass. (0->25)
	t17 = cos(qJ(4));
	t15 = sin(qJ(4));
	t30 = pkin(8) + r_i_i_C(2);
	t22 = t30 * t15;
	t34 = t17 * pkin(4) + pkin(3) + t22;
	t33 = pkin(1) + pkin(2);
	t14 = sin(qJ(5));
	t16 = cos(qJ(5));
	t25 = r_i_i_C(3) + qJ(6);
	t31 = pkin(5) + r_i_i_C(1);
	t32 = t25 * t14 + t31 * t16 + pkin(4);
	t29 = cos(qJ(1));
	t28 = sin(qJ(1));
	t27 = t14 * t17;
	t26 = t16 * t17;
	t24 = cos(pkin(9));
	t23 = sin(pkin(9));
	t10 = t29 * t23 - t28 * t24;
	t9 = -t28 * t23 - t29 * t24;
	t20 = t10 * t26 + t9 * t14;
	t19 = t10 * t27 - t9 * t16;
	t18 = t32 * t15 - t30 * t17;
	t6 = t10 * t14 - t9 * t26;
	t5 = -t10 * t16 - t9 * t27;
	t1 = [t9 * pkin(7) + t29 * qJ(2) + t34 * t10 + t25 * t19 + t31 * t20 - t33 * t28, t28, 0, t18 * t9, t25 * t6 - t31 * t5, t5; t10 * pkin(7) + t28 * qJ(2) + t25 * t5 + t33 * t29 + t31 * t6 - t34 * t9, -t29, 0, t18 * t10, t19 * t31 - t25 * t20, -t19; 0, 0, -1, -t32 * t17 - t22, (t31 * t14 - t25 * t16) * t15, -t15 * t14;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,6);
end