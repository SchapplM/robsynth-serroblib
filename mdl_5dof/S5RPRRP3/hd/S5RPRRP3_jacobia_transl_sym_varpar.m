% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRRP3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRP3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
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
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(8);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, 1, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(1)) * pkin(1), 0, 0, 0, 0; t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + sin(qJ(1)) * pkin(1), 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (27->10), mult. (26->12), div. (0->0), fcn. (26->6), ass. (0->10)
	t9 = pkin(6) + r_i_i_C(3);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t7 = r_i_i_C(1) * t4 + r_i_i_C(2) * t5;
	t6 = pkin(2) + t8;
	t3 = qJ(1) + pkin(8);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [0, 1, t8, 0, 0; cos(qJ(1)) * pkin(1) + t9 * t1 + t6 * t2, 0, -t7 * t1, 0, 0; sin(qJ(1)) * pkin(1) - t9 * t2 + t6 * t1, 0, t7 * t2, 0, 0;];
	Ja_transl = t10;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (59->14), mult. (43->16), div. (0->0), fcn. (43->8), ass. (0->14)
	t11 = qJ(3) + qJ(4);
	t7 = sin(t11);
	t8 = cos(t11);
	t24 = r_i_i_C(1) * t7 + r_i_i_C(2) * t8;
	t16 = t8 * r_i_i_C(1) - t7 * r_i_i_C(2);
	t23 = t16 + cos(qJ(3)) * pkin(3);
	t10 = qJ(1) + pkin(8);
	t6 = cos(t10);
	t19 = t24 * t6;
	t18 = sin(qJ(3)) * pkin(3);
	t17 = r_i_i_C(3) + pkin(7) + pkin(6);
	t14 = pkin(2) + t23;
	t5 = sin(t10);
	t1 = [0, 1, t23, t16, 0; cos(qJ(1)) * pkin(1) + t17 * t5 + t14 * t6, 0, (-t24 - t18) * t5, -t24 * t5, 0; sin(qJ(1)) * pkin(1) - t17 * t6 + t14 * t5, 0, t6 * t18 + t19, t19, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:48:13
	% EndTime: 2020-01-03 11:48:13
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (83->20), mult. (53->21), div. (0->0), fcn. (55->8), ass. (0->16)
	t23 = pkin(4) + r_i_i_C(1);
	t14 = qJ(3) + qJ(4);
	t10 = cos(t14);
	t9 = sin(t14);
	t16 = -t9 * r_i_i_C(2) + t23 * t10;
	t22 = cos(qJ(3)) * pkin(3) + t16;
	t13 = qJ(1) + pkin(8);
	t8 = cos(t13);
	t21 = t8 * t9;
	t18 = r_i_i_C(2) * t10;
	t19 = r_i_i_C(1) * t21 + t8 * t18;
	t17 = r_i_i_C(3) + qJ(5) + pkin(7) + pkin(6);
	t15 = pkin(2) + t22;
	t7 = sin(t13);
	t2 = -pkin(4) * t9 - sin(qJ(3)) * pkin(3);
	t1 = [0, 1, t22, t16, 0; cos(qJ(1)) * pkin(1) + t17 * t7 + t15 * t8, 0, (-r_i_i_C(1) * t9 - t18 + t2) * t7, (-t23 * t9 - t18) * t7, -t8; sin(qJ(1)) * pkin(1) - t17 * t8 + t15 * t7, 0, -t8 * t2 + t19, pkin(4) * t21 + t19, -t7;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end