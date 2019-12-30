% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPP1
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPPP1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPP1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPP1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:40
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:40
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:40
	% EndTime: 2019-12-29 18:08:40
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(7) + r_i_i_C(3);
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
	% StartTime: 2019-12-29 18:08:40
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (45->23), mult. (118->38), div. (0->0), fcn. (139->8), ass. (0->19)
	t6 = cos(pkin(5));
	t7 = sin(qJ(2));
	t20 = t6 * t7;
	t3 = sin(pkin(8));
	t5 = cos(pkin(8));
	t9 = cos(qJ(2));
	t12 = -(t3 * t20 - t5 * t9) * r_i_i_C(1) - (t5 * t20 + t3 * t9) * r_i_i_C(2) + t9 * pkin(2);
	t4 = sin(pkin(5));
	t21 = t4 * t7;
	t22 = qJ(3) * t21 + pkin(1) + t12;
	t19 = t6 * t9;
	t16 = t4 * (r_i_i_C(3) + qJ(3));
	t13 = t6 * qJ(3) + pkin(7) + (r_i_i_C(1) * t3 + r_i_i_C(2) * t5) * t4;
	t11 = (-t3 * t19 - t5 * t7) * r_i_i_C(1) + (-t5 * t19 + t3 * t7) * r_i_i_C(2) - t7 * pkin(2) + t9 * t16;
	t10 = cos(qJ(1));
	t8 = sin(qJ(1));
	t2 = t10 * t21 + t8 * t6;
	t1 = -t10 * t6 + t8 * t21;
	t14 = [-t1 * r_i_i_C(3) + t13 * t10 - t22 * t8, t11 * t10, t2, 0, 0; t2 * r_i_i_C(3) + t22 * t10 + t13 * t8, t11 * t8, t1, 0, 0; 0, t7 * t16 + t12, -t9 * t4, 0, 0;];
	Ja_transl = t14;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:41
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (74->33), mult. (199->50), div. (0->0), fcn. (242->8), ass. (0->32)
	t16 = cos(pkin(5));
	t17 = sin(qJ(2));
	t19 = cos(qJ(2));
	t15 = cos(pkin(8));
	t33 = t19 * t15;
	t13 = sin(pkin(8));
	t36 = t17 * t13;
	t27 = -t16 * t33 + t36;
	t14 = sin(pkin(5));
	t29 = t14 * (r_i_i_C(1) + qJ(3));
	t30 = r_i_i_C(3) + qJ(4);
	t34 = t19 * t13;
	t37 = t15 * t17;
	t41 = pkin(3) - r_i_i_C(2);
	t42 = -t30 * t27 - t41 * (t16 * t34 + t37) - pkin(2) * t17 + t19 * t29;
	t40 = t19 * pkin(2);
	t18 = sin(qJ(1));
	t39 = t14 * t18;
	t20 = cos(qJ(1));
	t38 = t14 * t20;
	t35 = t18 * t16;
	t32 = t19 * t20;
	t31 = t20 * t16;
	t28 = qJ(3) * t16 + pkin(7);
	t26 = t17 * t31 - t39;
	t25 = t17 * t35 + t38;
	t24 = t17 * t14 * qJ(3) + pkin(1) + t40;
	t12 = t17 * t38 + t35;
	t11 = t17 * t39 - t31;
	t3 = t13 * t32 + t26 * t15;
	t1 = t25 * t15 + t18 * t34;
	t2 = [-t11 * r_i_i_C(1) + t28 * t20 - t41 * (-t25 * t13 + t18 * t33) - t30 * t1 - t24 * t18, t42 * t20, t12, t3, 0; t12 * r_i_i_C(1) + t41 * (-t26 * t13 + t15 * t32) + t30 * t3 + t28 * t18 + t24 * t20, t42 * t18, t11, t1, 0; 0, t40 + t30 * (t16 * t37 + t34) + t17 * t29 + t41 * (-t16 * t36 + t33), -t19 * t14, t27, 0;];
	Ja_transl = t2;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:42
	% EndTime: 2019-12-29 18:08:42
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (98->34), mult. (264->50), div. (0->0), fcn. (324->8), ass. (0->36)
	t19 = cos(pkin(5));
	t18 = cos(pkin(8));
	t22 = cos(qJ(2));
	t38 = t22 * t18;
	t16 = sin(pkin(8));
	t20 = sin(qJ(2));
	t42 = t20 * t16;
	t48 = -t19 * t42 + t38;
	t39 = t22 * t16;
	t41 = t20 * t18;
	t29 = t19 * t39 + t41;
	t30 = -t19 * t38 + t42;
	t17 = sin(pkin(5));
	t46 = pkin(4) + r_i_i_C(1);
	t31 = t17 * (qJ(3) + t46);
	t34 = pkin(3) + r_i_i_C(3) + qJ(5);
	t35 = r_i_i_C(2) + qJ(4);
	t47 = -pkin(2) * t20 + t22 * t31 - t34 * t29 - t35 * t30;
	t45 = t22 * pkin(2);
	t21 = sin(qJ(1));
	t44 = t17 * t21;
	t23 = cos(qJ(1));
	t43 = t17 * t23;
	t40 = t21 * t19;
	t37 = t22 * t23;
	t36 = t23 * t19;
	t32 = qJ(3) * t19 + pkin(7);
	t28 = t20 * t36 - t44;
	t27 = t20 * t17 * qJ(3) + pkin(1) + t45;
	t12 = t20 * t43 + t40;
	t11 = t20 * t44 - t36;
	t4 = -t28 * t16 + t18 * t37;
	t3 = t16 * t37 + t28 * t18;
	t2 = -t16 * t43 + t48 * t21;
	t1 = t21 * t39 + (t20 * t40 + t43) * t18;
	t5 = [-t35 * t1 - t46 * t11 - t34 * t2 - t27 * t21 + t32 * t23, t47 * t23, t12, t3, t4; t46 * t12 + t32 * t21 + t27 * t23 + t35 * t3 + t34 * t4, t47 * t21, t11, t1, t2; 0, t45 + t35 * (t19 * t41 + t39) + t20 * t31 + t34 * t48, -t22 * t17, t30, t29;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end