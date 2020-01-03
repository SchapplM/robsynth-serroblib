% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR7
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
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPPR7_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR7_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR7_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
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
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
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
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (29->10), mult. (31->12), div. (0->0), fcn. (33->6), ass. (0->10)
	t5 = qJ(2) + pkin(8);
	t2 = sin(t5);
	t3 = cos(t5);
	t14 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2) + cos(qJ(2)) * pkin(2);
	t13 = r_i_i_C(3) + qJ(3) + pkin(6);
	t11 = pkin(1) + t14;
	t10 = -sin(qJ(2)) * pkin(2) - r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t11 * t8 + t13 * t9, t10 * t9, t8, 0, 0; t11 * t9 + t13 * t8, t10 * t8, -t9, 0, 0; 0, t14, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (53->13), mult. (51->14), div. (0->0), fcn. (56->6), ass. (0->12)
	t12 = r_i_i_C(3) + qJ(4);
	t14 = pkin(3) - r_i_i_C(2);
	t5 = qJ(2) + pkin(8);
	t2 = sin(t5);
	t3 = cos(t5);
	t16 = cos(qJ(2)) * pkin(2) + t12 * t2 + t14 * t3;
	t15 = pkin(1) + t16;
	t13 = r_i_i_C(1) + qJ(3) + pkin(6);
	t10 = -sin(qJ(2)) * pkin(2) + t12 * t3 - t14 * t2;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [t13 * t9 - t15 * t8, t10 * t9, t8, t9 * t2, 0; t13 * t8 + t15 * t9, t10 * t8, -t9, t8 * t2, 0; 0, t16, 0, -t3, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:37:15
	% EndTime: 2019-12-31 19:37:15
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (90->26), mult. (103->36), div. (0->0), fcn. (116->8), ass. (0->22)
	t20 = pkin(3) + pkin(7) + r_i_i_C(3);
	t9 = qJ(2) + pkin(8);
	t7 = cos(t9);
	t27 = t20 * t7 + cos(qJ(2)) * pkin(2);
	t6 = sin(t9);
	t26 = t6 * qJ(4) + pkin(1) + t27;
	t25 = pkin(4) + qJ(3) + pkin(6);
	t11 = sin(qJ(5));
	t15 = cos(qJ(1));
	t24 = t11 * t15;
	t13 = sin(qJ(1));
	t23 = t13 * t11;
	t14 = cos(qJ(5));
	t22 = t13 * t14;
	t21 = t14 * t15;
	t17 = r_i_i_C(1) * t11 + r_i_i_C(2) * t14 + qJ(4);
	t16 = -sin(qJ(2)) * pkin(2) + t17 * t7 - t20 * t6;
	t4 = -t6 * t23 + t21;
	t3 = t6 * t22 + t24;
	t2 = t6 * t24 + t22;
	t1 = t6 * t21 - t23;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t26 * t13 + t25 * t15, t16 * t15, t13, t15 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t25 * t13 + t26 * t15, t16 * t13, -t15, t13 * t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, t17 * t6 + t27, 0, -t7, (-r_i_i_C(1) * t14 + r_i_i_C(2) * t11) * t7;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end