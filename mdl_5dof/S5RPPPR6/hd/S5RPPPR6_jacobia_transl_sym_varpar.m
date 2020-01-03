% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPPPR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPPR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
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
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(7)) - r_i_i_C(2) * sin(pkin(7)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t3 * t5 + t4 * t6, t3, 0, 0, 0; t3 * t6 + t4 * t5, -t4, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (14->9), mult. (24->8), div. (0->0), fcn. (29->4), ass. (0->7)
	t1 = sin(pkin(7));
	t2 = cos(pkin(7));
	t8 = pkin(1) + (pkin(2) - r_i_i_C(2)) * t2 + (r_i_i_C(3) + qJ(3)) * t1;
	t6 = r_i_i_C(1) + qJ(2);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t5 = [-t8 * t3 + t6 * t4, t3, t4 * t1, 0, 0; t6 * t3 + t8 * t4, -t4, t3 * t1, 0, 0; 0, 0, -t2, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (22->13), mult. (42->14), div. (0->0), fcn. (54->6), ass. (0->9)
	t1 = sin(pkin(8));
	t2 = sin(pkin(7));
	t3 = cos(pkin(8));
	t4 = cos(pkin(7));
	t10 = pkin(1) + (pkin(2) + r_i_i_C(3) + qJ(4)) * t4 + (r_i_i_C(1) * t1 + r_i_i_C(2) * t3 + qJ(3)) * t2;
	t7 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2) + pkin(3) + qJ(2);
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t8 = [-t10 * t5 + t7 * t6, t5, t6 * t2, t6 * t4, 0; t10 * t6 + t7 * t5, -t6, t5 * t2, t5 * t4, 0; 0, 0, -t4, t2, 0;];
	Ja_transl = t8;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:48:16
	% EndTime: 2019-12-31 17:48:16
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (49->33), mult. (110->48), div. (0->0), fcn. (142->8), ass. (0->24)
	t26 = r_i_i_C(3) + pkin(6);
	t11 = cos(pkin(7));
	t8 = sin(pkin(8));
	t25 = t11 * t8;
	t13 = sin(qJ(1));
	t9 = sin(pkin(7));
	t24 = t13 * t9;
	t15 = cos(qJ(1));
	t23 = t15 * t9;
	t10 = cos(pkin(8));
	t22 = t10 * t15;
	t21 = t13 * t10;
	t20 = t13 * t11;
	t19 = t15 * t11;
	t18 = pkin(2) + qJ(4);
	t17 = pkin(3) + qJ(2);
	t16 = -qJ(3) * t9 - pkin(1);
	t14 = cos(qJ(5));
	t12 = sin(qJ(5));
	t6 = -t8 * t24 + t22;
	t4 = t8 * t23 + t21;
	t2 = t12 * t19 + t4 * t14;
	t1 = -t4 * t12 + t14 * t19;
	t3 = [t26 * (t15 * t8 + t9 * t21) + t17 * t15 + (t14 * r_i_i_C(1) - t12 * r_i_i_C(2) + pkin(4)) * t6 + ((-r_i_i_C(1) * t12 - r_i_i_C(2) * t14 - t18) * t11 + t16) * t13, t13, t23, t19, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t4 * pkin(4) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t26 * (t13 * t8 - t9 * t22) + t17 * t13 + (t18 * t11 - t16) * t15, -t15, t24, t20, (t12 * t6 + t14 * t20) * r_i_i_C(1) + (-t12 * t20 + t14 * t6) * r_i_i_C(2); 0, 0, -t11, t9, (t12 * t25 + t14 * t9) * r_i_i_C(1) + (-t12 * t9 + t14 * t25) * r_i_i_C(2);];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,5);
end