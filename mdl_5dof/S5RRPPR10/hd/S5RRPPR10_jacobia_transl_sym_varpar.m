% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPPR10_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR10_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR10_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:31
	% EndTime: 2019-12-31 19:45:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:31
	% EndTime: 2019-12-31 19:45:32
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
	% StartTime: 2019-12-31 19:45:31
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(6) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t2 * t5 + t4 * t8, t6 * t4, 0, 0, 0; t2 * t8 + t4 * t5, t6 * t2, 0, 0, 0; 0, t7, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:31
	% EndTime: 2019-12-31 19:45:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (29->11), mult. (65->16), div. (0->0), fcn. (72->6), ass. (0->13)
	t1 = sin(pkin(8));
	t2 = cos(pkin(8));
	t10 = r_i_i_C(1) * t2 - r_i_i_C(2) * t1 + pkin(2);
	t11 = r_i_i_C(3) + qJ(3);
	t3 = sin(qJ(2));
	t5 = cos(qJ(2));
	t8 = t10 * t5 + t11 * t3;
	t12 = pkin(1) + t8;
	t9 = t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + pkin(6);
	t7 = -t10 * t3 + t11 * t5;
	t6 = cos(qJ(1));
	t4 = sin(qJ(1));
	t13 = [-t12 * t4 + t9 * t6, t7 * t6, t6 * t3, 0, 0; t12 * t6 + t9 * t4, t7 * t4, t4 * t3, 0, 0; 0, t8, -t5, 0, 0;];
	Ja_transl = t13;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:31
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (45->20), mult. (104->30), div. (0->0), fcn. (120->6), ass. (0->18)
	t15 = r_i_i_C(2) + qJ(3);
	t7 = sin(qJ(2));
	t12 = t15 * t7;
	t9 = cos(qJ(2));
	t20 = t9 * pkin(2) + pkin(1) + t12;
	t14 = r_i_i_C(3) + qJ(4);
	t17 = pkin(3) + r_i_i_C(1);
	t5 = sin(pkin(8));
	t6 = cos(pkin(8));
	t19 = t14 * t5 + t17 * t6 + pkin(2);
	t8 = sin(qJ(1));
	t18 = t8 * t9;
	t10 = cos(qJ(1));
	t16 = t10 * t9;
	t11 = t15 * t9 - t19 * t7;
	t3 = t5 * t16 - t8 * t6;
	t1 = t10 * t6 + t5 * t18;
	t2 = [pkin(6) * t10 + t17 * (t10 * t5 - t6 * t18) - t14 * t1 - t20 * t8, t11 * t10, t10 * t7, t3, 0; t8 * pkin(6) + t17 * (t6 * t16 + t8 * t5) + t14 * t3 + t20 * t10, t11 * t8, t8 * t7, t1, 0; 0, t19 * t9 + t12, -t9, t7 * t5, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:45:31
	% EndTime: 2019-12-31 19:45:32
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (88->34), mult. (217->55), div. (0->0), fcn. (263->8), ass. (0->25)
	t14 = cos(qJ(2));
	t11 = sin(qJ(2));
	t21 = -r_i_i_C(3) - pkin(7) + qJ(3);
	t19 = t21 * t11;
	t26 = t14 * pkin(2) + pkin(1) + t19;
	t10 = sin(qJ(5));
	t13 = cos(qJ(5));
	t24 = pkin(3) + pkin(4);
	t17 = t13 * r_i_i_C(1) - t10 * r_i_i_C(2) + t24;
	t18 = t10 * r_i_i_C(1) + t13 * r_i_i_C(2) + qJ(4);
	t8 = sin(pkin(8));
	t9 = cos(pkin(8));
	t25 = t17 * t9 + t18 * t8 + pkin(2);
	t12 = sin(qJ(1));
	t23 = t12 * t14;
	t15 = cos(qJ(1));
	t22 = t14 * t15;
	t16 = -t25 * t11 + t21 * t14;
	t6 = t12 * t8 + t9 * t22;
	t5 = -t12 * t9 + t8 * t22;
	t4 = -t15 * t8 + t9 * t23;
	t3 = t15 * t9 + t8 * t23;
	t2 = t10 * t5 + t13 * t6;
	t1 = -t10 * t6 + t13 * t5;
	t7 = [t15 * pkin(6) - t26 * t12 - t17 * t4 - t18 * t3, t16 * t15, t15 * t11, t5, r_i_i_C(1) * t1 - t2 * r_i_i_C(2); t12 * pkin(6) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * qJ(4) + t26 * t15 + t24 * t6, t16 * t12, t12 * t11, t3, (-t10 * t4 + t13 * t3) * r_i_i_C(1) + (-t10 * t3 - t13 * t4) * r_i_i_C(2); 0, t25 * t14 + t19, -t14, t11 * t8, ((-t10 * t9 + t13 * t8) * r_i_i_C(1) + (-t10 * t8 - t13 * t9) * r_i_i_C(2)) * t11;];
	Ja_transl = t7;
else
	Ja_transl=NaN(3,5);
end