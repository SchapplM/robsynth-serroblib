% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPPRR8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRR8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:01:29
	% EndTime: 2019-12-31 18:01:30
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:01:29
	% EndTime: 2019-12-31 18:01:30
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
	% StartTime: 2019-12-31 18:01:29
	% EndTime: 2019-12-31 18:01:30
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (7->5), mult. (8->4), div. (0->0), fcn. (10->2), ass. (0->5)
	t4 = pkin(1) + r_i_i_C(1);
	t3 = r_i_i_C(3) + qJ(2);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t4 * t1 + t3 * t2, t1, 0, 0, 0; t3 * t1 + t4 * t2, -t2, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:01:29
	% EndTime: 2019-12-31 18:01:30
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (13->10), mult. (18->12), div. (0->0), fcn. (24->4), ass. (0->8)
	t7 = pkin(1) + pkin(2);
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t4 = cos(pkin(8));
	t3 = sin(pkin(8));
	t2 = t5 * t3 + t6 * t4;
	t1 = t6 * t3 - t5 * t4;
	t8 = [t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + t6 * qJ(2) - t7 * t5, t5, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + t5 * qJ(2) + t7 * t6, -t6, 0, 0, 0; 0, 0, -1, 0, 0;];
	Ja_transl = t8;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:01:29
	% EndTime: 2019-12-31 18:01:29
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (39->13), mult. (36->14), div. (0->0), fcn. (46->6), ass. (0->12)
	t16 = pkin(1) + cos(pkin(8)) * pkin(3) + pkin(2);
	t15 = pkin(8) + qJ(4);
	t14 = pkin(3) * sin(pkin(8)) + qJ(2);
	t13 = cos(t15);
	t12 = sin(t15);
	t8 = sin(qJ(1));
	t9 = cos(qJ(1));
	t1 = -t8 * t12 - t9 * t13;
	t2 = t9 * t12 - t8 * t13;
	t11 = -t2 * r_i_i_C(1) + t1 * r_i_i_C(2);
	t10 = t1 * r_i_i_C(1) + t2 * r_i_i_C(2);
	t3 = [t14 * t9 - t16 * t8 - t11, t8, 0, t11, 0; t14 * t8 + t16 * t9 - t10, -t9, 0, t10, 0; 0, 0, -1, 0, 0;];
	Ja_transl = t3;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:01:29
	% EndTime: 2019-12-31 18:01:30
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (102->17), mult. (98->20), div. (0->0), fcn. (124->8), ass. (0->18)
	t23 = pkin(1) + cos(pkin(8)) * pkin(3) + pkin(2);
	t22 = qJ(2) + pkin(3) * sin(pkin(8));
	t15 = pkin(8) + qJ(4);
	t13 = sin(t15);
	t14 = cos(t15);
	t16 = sin(qJ(1));
	t17 = cos(qJ(1));
	t1 = -t16 * t13 - t17 * t14;
	t8 = sin(qJ(5));
	t9 = cos(qJ(5));
	t12 = -r_i_i_C(1) * t9 + t8 * r_i_i_C(2);
	t10 = pkin(4) - t12;
	t18 = pkin(7) + r_i_i_C(3);
	t2 = t17 * t13 - t16 * t14;
	t21 = -t18 * t1 - t10 * t2;
	t20 = t10 * t1 - t18 * t2;
	t11 = r_i_i_C(1) * t8 + r_i_i_C(2) * t9;
	t3 = [-t23 * t16 + t22 * t17 - t21, t16, 0, t21, t11 * t1; t22 * t16 + t23 * t17 - t20, -t17, 0, t20, t11 * t2; 0, 0, -1, 0, t12;];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,5);
end