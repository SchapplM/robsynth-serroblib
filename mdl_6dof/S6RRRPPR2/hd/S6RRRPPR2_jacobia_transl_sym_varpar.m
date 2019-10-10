% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRRPPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% Ja_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobia_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobia_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
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
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(7) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t2 * t5 + t4 * t8, t6 * t4, 0, 0, 0, 0; t2 * t8 + t4 * t5, t6 * t2, 0, 0, 0, 0; 0, t7, 0, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (37->9), mult. (41->14), div. (0->0), fcn. (41->6), ass. (0->12)
	t6 = qJ(2) + qJ(3);
	t3 = sin(t6);
	t4 = cos(t6);
	t14 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t17 = t14 + cos(qJ(2)) * pkin(2);
	t15 = r_i_i_C(3) + pkin(8) + pkin(7);
	t13 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t12 = pkin(1) + t17;
	t11 = -sin(qJ(2)) * pkin(2) + t13;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t12 * t8 + t15 * t9, t11 * t9, t13 * t9, 0, 0, 0; t12 * t9 + t15 * t8, t11 * t8, t13 * t8, 0, 0, 0; 0, t17, t14, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (72->13), mult. (51->16), div. (0->0), fcn. (53->8), ass. (0->13)
	t11 = qJ(2) + qJ(3);
	t7 = pkin(10) + t11;
	t4 = sin(t7);
	t5 = cos(t7);
	t18 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2) + pkin(3) * cos(t11);
	t22 = t18 + cos(qJ(2)) * pkin(2);
	t14 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5 - pkin(3) * sin(t11);
	t19 = r_i_i_C(3) + qJ(4) + pkin(8) + pkin(7);
	t16 = pkin(1) + t22;
	t15 = -sin(qJ(2)) * pkin(2) + t14;
	t13 = cos(qJ(1));
	t12 = sin(qJ(1));
	t1 = [-t12 * t16 + t13 * t19, t15 * t13, t14 * t13, t12, 0, 0; t12 * t19 + t13 * t16, t15 * t12, t14 * t12, -t13, 0, 0; 0, t22, t18, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (127->21), mult. (81->22), div. (0->0), fcn. (86->8), ass. (0->19)
	t39 = r_i_i_C(3) + qJ(5);
	t19 = qJ(2) + qJ(3);
	t15 = pkin(10) + t19;
	t12 = sin(t15);
	t13 = cos(t15);
	t22 = pkin(3) * cos(t19) + (pkin(4) - r_i_i_C(2)) * t13 + t39 * t12;
	t38 = cos(qJ(2)) * pkin(2) + t22;
	t37 = t39 * t13;
	t23 = -pkin(4) * t12 - pkin(3) * sin(t19);
	t36 = pkin(1) + t38;
	t31 = r_i_i_C(1) + qJ(4) + pkin(8) + pkin(7);
	t20 = sin(qJ(1));
	t30 = t20 * t12;
	t21 = cos(qJ(1));
	t29 = t21 * t12;
	t26 = r_i_i_C(2) * t30 + t37 * t20;
	t25 = r_i_i_C(2) * t29 + t37 * t21;
	t24 = -sin(qJ(2)) * pkin(2) + t23;
	t1 = [-t36 * t20 + t31 * t21, t24 * t21 + t25, t23 * t21 + t25, t20, t29, 0; t31 * t20 + t36 * t21, t24 * t20 + t26, t23 * t20 + t26, -t21, t30, 0; 0, t38, t22, 0, -t13, 0;];
	Ja_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobia_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (198->32), mult. (149->42), div. (0->0), fcn. (162->10), ass. (0->29)
	t27 = sin(qJ(6));
	t29 = cos(qJ(6));
	t55 = r_i_i_C(1) * t27 + r_i_i_C(2) * t29;
	t26 = qJ(2) + qJ(3);
	t22 = pkin(10) + t26;
	t19 = sin(t22);
	t20 = cos(t22);
	t40 = pkin(4) + pkin(9) + r_i_i_C(3);
	t54 = t40 * t20 + t19 * qJ(5) + pkin(3) * cos(t26);
	t52 = (qJ(5) + t55) * t20;
	t31 = -t40 * t19 - pkin(3) * sin(t26);
	t24 = cos(qJ(2)) * pkin(2);
	t50 = pkin(1) + t24 + t54;
	t46 = pkin(5) + qJ(4) + pkin(8) + pkin(7);
	t28 = sin(qJ(1));
	t45 = t28 * t27;
	t44 = t28 * t29;
	t30 = cos(qJ(1));
	t43 = t30 * t27;
	t42 = t30 * t29;
	t39 = t52 * t28;
	t38 = t52 * t30;
	t33 = -sin(qJ(2)) * pkin(2) + t31;
	t32 = t55 * t19 + t54;
	t4 = -t19 * t45 + t42;
	t3 = t19 * t44 + t43;
	t2 = t19 * t43 + t44;
	t1 = t19 * t42 - t45;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t50 * t28 + t46 * t30, t33 * t30 + t38, t31 * t30 + t38, t28, t30 * t19, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t46 * t28 + t50 * t30, t33 * t28 + t39, t31 * t28 + t39, -t30, t28 * t19, t3 * r_i_i_C(1) + t4 * r_i_i_C(2); 0, t24 + t32, t32, 0, -t20, (-r_i_i_C(1) * t29 + r_i_i_C(2) * t27) * t20;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,6);
end