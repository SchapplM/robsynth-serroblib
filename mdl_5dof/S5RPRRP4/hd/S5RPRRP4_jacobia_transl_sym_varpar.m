% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRP4
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRRP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
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
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (10->7), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t3 * t6 + t4 * t5, -t4, 0, 0, 0; t3 * t5 - t4 * t6, -t3, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (25->18), mult. (52->25), div. (0->0), fcn. (62->6), ass. (0->16)
	t7 = sin(qJ(3));
	t8 = sin(qJ(1));
	t15 = t8 * t7;
	t9 = cos(qJ(3));
	t14 = t8 * t9;
	t10 = cos(qJ(1));
	t13 = t10 * t7;
	t12 = t10 * t9;
	t5 = sin(pkin(8));
	t6 = cos(pkin(8));
	t11 = pkin(2) * t6 + pkin(1) + (pkin(6) + r_i_i_C(3)) * t5;
	t4 = t6 * t12 + t15;
	t3 = t6 * t13 - t14;
	t2 = t6 * t14 - t13;
	t1 = -t6 * t15 - t12;
	t16 = [0, 0, (-r_i_i_C(1) * t7 - r_i_i_C(2) * t9) * t5, 0, 0; t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t8 * qJ(2) + t11 * t10, -t10, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t10 * qJ(2) + t11 * t8, -t8, t3 * r_i_i_C(1) + t4 * r_i_i_C(2), 0, 0;];
	Ja_transl = t16;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (71->27), mult. (90->36), div. (0->0), fcn. (104->8), ass. (0->22)
	t12 = qJ(3) + qJ(4);
	t10 = sin(t12);
	t11 = cos(t12);
	t18 = cos(qJ(1));
	t14 = cos(pkin(8));
	t16 = sin(qJ(1));
	t24 = t14 * t16;
	t5 = -t10 * t24 - t11 * t18;
	t6 = -t10 * t18 + t11 * t24;
	t27 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t23 = t14 * t18;
	t7 = t10 * t23 - t11 * t16;
	t8 = t10 * t16 + t11 * t23;
	t26 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t15 = sin(qJ(3));
	t25 = pkin(3) * t15;
	t22 = qJ(2) + t25;
	t21 = -r_i_i_C(1) * t10 - r_i_i_C(2) * t11;
	t13 = sin(pkin(8));
	t17 = cos(qJ(3));
	t20 = t14 * (pkin(3) * t17 + pkin(2)) + pkin(1) + (r_i_i_C(3) + pkin(7) + pkin(6)) * t13;
	t1 = [0, 0, (t21 - t25) * t13, t21 * t13, 0; t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t22 * t16 + t20 * t18, -t18, (-t15 * t24 - t17 * t18) * pkin(3) + t27, t27, 0; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t20 * t16 - t22 * t18, -t16, (t15 * t23 - t16 * t17) * pkin(3) + t26, t26, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (103->36), mult. (113->41), div. (0->0), fcn. (130->8), ass. (0->22)
	t17 = qJ(3) + qJ(4);
	t14 = cos(t17);
	t19 = cos(pkin(8));
	t21 = cos(qJ(1));
	t13 = sin(t17);
	t20 = sin(qJ(1));
	t25 = t20 * t13;
	t5 = -t14 * t21 - t19 * t25;
	t24 = t20 * t14;
	t6 = -t13 * t21 + t19 * t24;
	t29 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t26 = t19 * t21;
	t7 = t13 * t26 - t24;
	t8 = t14 * t26 + t25;
	t28 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
	t27 = r_i_i_C(2) * t14;
	t11 = pkin(4) * t14 + cos(qJ(3)) * pkin(3);
	t10 = pkin(4) * t13 + sin(qJ(3)) * pkin(3);
	t23 = qJ(2) + t10;
	t18 = sin(pkin(8));
	t22 = t19 * (pkin(2) + t11) + pkin(1) + (r_i_i_C(3) + qJ(5) + pkin(7) + pkin(6)) * t18;
	t1 = [0, 0, (-r_i_i_C(1) * t13 - t10 - t27) * t18, (-t27 + (-pkin(4) - r_i_i_C(1)) * t13) * t18, -t19; t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t23 * t20 + t22 * t21, -t21, -t20 * t19 * t10 - t11 * t21 + t29, t5 * pkin(4) + t29, t20 * t18; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t22 * t20 - t23 * t21, -t20, t10 * t26 - t20 * t11 + t28, t7 * pkin(4) + t28, -t21 * t18;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end