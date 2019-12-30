% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:30
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRRP5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP5_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRP5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP5_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:30
	% EndTime: 2019-12-29 20:30:30
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:25
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:25
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(6) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t5 * t2 + t8 * t4, t6 * t4, 0, 0, 0; t8 * t2 + t5 * t4, t6 * t2, 0, 0, 0; 0, t7, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:25
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (37->9), mult. (41->14), div. (0->0), fcn. (41->6), ass. (0->12)
	t6 = qJ(2) + qJ(3);
	t3 = sin(t6);
	t4 = cos(t6);
	t14 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t17 = t14 + cos(qJ(2)) * pkin(2);
	t15 = r_i_i_C(3) + pkin(7) + pkin(6);
	t13 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t12 = pkin(1) + t17;
	t11 = -sin(qJ(2)) * pkin(2) + t13;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t12 * t8 + t15 * t9, t11 * t9, t13 * t9, 0, 0; t12 * t9 + t15 * t8, t11 * t8, t13 * t8, 0, 0; 0, t17, t14, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:25
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (86->12), mult. (61->18), div. (0->0), fcn. (61->8), ass. (0->15)
	t11 = qJ(2) + qJ(3);
	t8 = qJ(4) + t11;
	t5 = sin(t8);
	t6 = cos(t8);
	t19 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t18 = t19 + pkin(3) * cos(t11);
	t23 = t18 + cos(qJ(2)) * pkin(2);
	t17 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t14 = t17 - pkin(3) * sin(t11);
	t20 = r_i_i_C(3) + pkin(8) + pkin(7) + pkin(6);
	t16 = pkin(1) + t23;
	t15 = -sin(qJ(2)) * pkin(2) + t14;
	t13 = cos(qJ(1));
	t12 = sin(qJ(1));
	t1 = [-t16 * t12 + t20 * t13, t15 * t13, t14 * t13, t17 * t13, 0; t20 * t12 + t16 * t13, t15 * t12, t14 * t12, t17 * t12, 0; 0, t23, t18, t19, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:30:25
	% EndTime: 2019-12-29 20:30:26
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (159->20), mult. (101->22), div. (0->0), fcn. (104->8), ass. (0->20)
	t38 = r_i_i_C(3) + qJ(5);
	t18 = qJ(2) + qJ(3);
	t15 = qJ(4) + t18;
	t12 = sin(t15);
	t13 = cos(t15);
	t33 = pkin(4) + r_i_i_C(1);
	t24 = t38 * t12 + t33 * t13;
	t23 = pkin(3) * cos(t18) + t24;
	t37 = cos(qJ(2)) * pkin(2) + t23;
	t36 = t38 * t13;
	t25 = t33 * t12;
	t21 = -t25 - pkin(3) * sin(t18);
	t34 = pkin(1) + t37;
	t19 = sin(qJ(1));
	t32 = t36 * t19;
	t20 = cos(qJ(1));
	t31 = t36 * t20;
	t28 = r_i_i_C(2) + pkin(8) + pkin(7) + pkin(6);
	t22 = -sin(qJ(2)) * pkin(2) + t21;
	t1 = [-t34 * t19 + t28 * t20, t22 * t20 + t31, t21 * t20 + t31, -t20 * t25 + t31, t20 * t12; t28 * t19 + t34 * t20, t22 * t19 + t32, t21 * t19 + t32, -t19 * t25 + t32, t19 * t12; 0, t37, t23, t24, -t13;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end