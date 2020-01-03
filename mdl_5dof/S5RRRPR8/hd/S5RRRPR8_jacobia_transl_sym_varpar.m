% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRPR8_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR8_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR8_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:21:20
	% EndTime: 2019-12-31 21:21:20
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:21:20
	% EndTime: 2019-12-31 21:21:20
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
	% StartTime: 2019-12-31 21:21:20
	% EndTime: 2019-12-31 21:21:20
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
	% StartTime: 2019-12-31 21:21:20
	% EndTime: 2019-12-31 21:21:20
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-12-31 21:21:20
	% EndTime: 2019-12-31 21:21:21
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (73->17), mult. (71->20), div. (0->0), fcn. (74->6), ass. (0->17)
	t33 = r_i_i_C(3) + qJ(4);
	t14 = qJ(2) + qJ(3);
	t11 = sin(t14);
	t12 = cos(t14);
	t19 = (pkin(3) - r_i_i_C(2)) * t12 + t33 * t11;
	t32 = cos(qJ(2)) * pkin(2) + t19;
	t31 = t33 * t12;
	t30 = pkin(1) + t32;
	t27 = r_i_i_C(1) + pkin(7) + pkin(6);
	t16 = sin(qJ(1));
	t26 = t16 * t11;
	t17 = cos(qJ(1));
	t25 = t17 * t11;
	t22 = r_i_i_C(2) * t26 + t31 * t16;
	t21 = r_i_i_C(2) * t25 + t31 * t17;
	t20 = -sin(qJ(2)) * pkin(2) - pkin(3) * t11;
	t1 = [-t30 * t16 + t27 * t17, t20 * t17 + t21, -pkin(3) * t25 + t21, t25, 0; t27 * t16 + t30 * t17, t20 * t16 + t22, -pkin(3) * t26 + t22, t26, 0; 0, t32, t19, -t12, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:21:20
	% EndTime: 2019-12-31 21:21:21
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (122->28), mult. (139->40), div. (0->0), fcn. (150->8), ass. (0->28)
	t22 = sin(qJ(5));
	t25 = cos(qJ(5));
	t49 = r_i_i_C(1) * t22 + r_i_i_C(2) * t25;
	t21 = qJ(2) + qJ(3);
	t18 = sin(t21);
	t19 = cos(t21);
	t36 = pkin(3) + pkin(8) + r_i_i_C(3);
	t48 = t18 * qJ(4) + t36 * t19;
	t46 = (qJ(4) + t49) * t19;
	t20 = cos(qJ(2)) * pkin(2);
	t45 = t20 + pkin(1) + t48;
	t42 = pkin(4) + pkin(7) + pkin(6);
	t24 = sin(qJ(1));
	t41 = t24 * t22;
	t40 = t24 * t25;
	t26 = cos(qJ(1));
	t39 = t26 * t22;
	t38 = t26 * t25;
	t35 = t46 * t24;
	t34 = t46 * t26;
	t30 = t36 * t18;
	t29 = t49 * t18 + t48;
	t28 = -sin(qJ(2)) * pkin(2) - t30;
	t4 = -t18 * t41 + t38;
	t3 = t18 * t40 + t39;
	t2 = t18 * t39 + t40;
	t1 = t18 * t38 - t41;
	t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t45 * t24 + t42 * t26, t28 * t26 + t34, -t26 * t30 + t34, t26 * t18, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t42 * t24 + t45 * t26, t28 * t24 + t35, -t24 * t30 + t35, t24 * t18, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, t20 + t29, t29, -t19, (-r_i_i_C(1) * t25 + r_i_i_C(2) * t22) * t19;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end