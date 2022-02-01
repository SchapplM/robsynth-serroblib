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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
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
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (24->17), mult. (48->25), div. (0->0), fcn. (58->6), ass. (0->16)
	t8 = sin(qJ(3));
	t9 = sin(qJ(1));
	t16 = t9 * t8;
	t11 = cos(qJ(1));
	t15 = t11 * t8;
	t10 = cos(qJ(3));
	t14 = t9 * t10;
	t13 = t11 * t10;
	t6 = sin(pkin(8));
	t7 = cos(pkin(8));
	t12 = pkin(2) * t7 + pkin(1) + (pkin(6) + r_i_i_C(3)) * t6;
	t4 = t7 * t13 + t16;
	t3 = -t7 * t15 + t14;
	t2 = -t7 * t14 + t15;
	t1 = t7 * t16 + t13;
	t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t11 * qJ(2) - t12 * t9, t9, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t9 * qJ(2) + t12 * t11, -t11, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0; 0, 0, (-r_i_i_C(1) * t8 - r_i_i_C(2) * t10) * t6, 0, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (70->26), mult. (86->36), div. (0->0), fcn. (100->8), ass. (0->23)
	t13 = qJ(3) + qJ(4);
	t12 = cos(t13);
	t15 = cos(pkin(8));
	t19 = cos(qJ(1));
	t11 = sin(t13);
	t17 = sin(qJ(1));
	t23 = t17 * t11;
	t6 = t12 * t19 + t15 * t23;
	t22 = t17 * t12;
	t7 = t11 * t19 - t15 * t22;
	t27 = -t6 * r_i_i_C(1) + t7 * r_i_i_C(2);
	t24 = t15 * t19;
	t8 = -t11 * t24 + t22;
	t9 = t12 * t24 + t23;
	t26 = t8 * r_i_i_C(1) - t9 * r_i_i_C(2);
	t16 = sin(qJ(3));
	t25 = pkin(3) * t16;
	t14 = sin(pkin(8));
	t18 = cos(qJ(3));
	t21 = pkin(1) + (pkin(3) * t18 + pkin(2)) * t15 + (r_i_i_C(3) + pkin(7) + pkin(6)) * t14;
	t20 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
	t10 = qJ(2) + t25;
	t1 = [t7 * r_i_i_C(1) + t6 * r_i_i_C(2) + t10 * t19 - t21 * t17, t17, (-t16 * t24 + t17 * t18) * pkin(3) + t26, t26, 0; t9 * r_i_i_C(1) + t8 * r_i_i_C(2) + t10 * t17 + t21 * t19, -t19, (-t15 * t16 * t17 - t18 * t19) * pkin(3) + t27, t27, 0; 0, 0, (t20 - t25) * t14, t20 * t14, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (101->34), mult. (113->41), div. (0->0), fcn. (130->8), ass. (0->22)
	t17 = qJ(3) + qJ(4);
	t14 = cos(t17);
	t11 = pkin(4) * t14 + cos(qJ(3)) * pkin(3);
	t18 = sin(pkin(8));
	t19 = cos(pkin(8));
	t31 = (r_i_i_C(3) + qJ(5) + pkin(7) + pkin(6)) * t18 + t19 * (pkin(2) + t11) + pkin(1);
	t21 = cos(qJ(1));
	t13 = sin(t17);
	t20 = sin(qJ(1));
	t25 = t20 * t13;
	t5 = t14 * t21 + t19 * t25;
	t24 = t20 * t14;
	t6 = t13 * t21 - t19 * t24;
	t30 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t26 = t19 * t21;
	t7 = -t13 * t26 + t24;
	t8 = t14 * t26 + t25;
	t29 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t28 = r_i_i_C(2) * t14;
	t10 = pkin(4) * t13 + sin(qJ(3)) * pkin(3);
	t23 = qJ(2) + t10;
	t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t31 * t20 + t23 * t21, t20, -t10 * t26 + t20 * t11 + t29, t7 * pkin(4) + t29, t21 * t18; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t23 * t20 + t31 * t21, -t21, -t20 * t19 * t10 - t11 * t21 + t30, -t5 * pkin(4) + t30, t20 * t18; 0, 0, (-r_i_i_C(1) * t13 - t10 - t28) * t18, (-t28 + (-pkin(4) - r_i_i_C(1)) * t13) * t18, -t19;];
	Ja_transl = t1;
end