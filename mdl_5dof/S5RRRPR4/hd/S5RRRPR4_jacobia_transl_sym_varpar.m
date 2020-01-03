% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR4
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRPR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:12:10
	% EndTime: 2019-12-31 21:12:10
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:12:10
	% EndTime: 2019-12-31 21:12:10
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
	% StartTime: 2019-12-31 21:12:10
	% EndTime: 2019-12-31 21:12:10
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t4 = qJ(1) + qJ(2);
	t2 = sin(t4);
	t3 = cos(t4);
	t6 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t5 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t1 = [-sin(qJ(1)) * pkin(1) + t5, t5, 0, 0, 0; cos(qJ(1)) * pkin(1) + t6, t6, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:12:10
	% EndTime: 2019-12-31 21:12:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (45->9), mult. (40->12), div. (0->0), fcn. (40->6), ass. (0->12)
	t20 = r_i_i_C(3) + pkin(7);
	t11 = sin(qJ(3));
	t12 = cos(qJ(3));
	t19 = t12 * r_i_i_C(1) - t11 * r_i_i_C(2);
	t18 = -pkin(2) - t19;
	t15 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
	t10 = qJ(1) + qJ(2);
	t8 = sin(t10);
	t9 = cos(t10);
	t14 = -t18 * t9 + t20 * t8;
	t13 = t18 * t8 + t20 * t9;
	t1 = [-sin(qJ(1)) * pkin(1) + t13, t13, t15 * t9, 0, 0; cos(qJ(1)) * pkin(1) + t14, t14, t15 * t8, 0, 0; 0, 0, t19, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:12:10
	% EndTime: 2019-12-31 21:12:10
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (74->14), mult. (68->15), div. (0->0), fcn. (71->6), ass. (0->15)
	t14 = cos(qJ(3));
	t21 = pkin(3) + r_i_i_C(1);
	t24 = t21 * t14;
	t23 = r_i_i_C(2) + pkin(7);
	t18 = r_i_i_C(3) + qJ(4);
	t13 = sin(qJ(3));
	t22 = t18 * t13 + t24;
	t12 = qJ(1) + qJ(2);
	t11 = cos(t12);
	t20 = t11 * t13;
	t10 = sin(t12);
	t17 = t23 * t10 + t18 * t20 + (pkin(2) + t24) * t11;
	t16 = -t21 * t13 + t18 * t14;
	t15 = (-pkin(2) - t22) * t10 + t23 * t11;
	t1 = [-sin(qJ(1)) * pkin(1) + t15, t15, t16 * t11, t20, 0; cos(qJ(1)) * pkin(1) + t17, t17, t16 * t10, t10 * t13, 0; 0, 0, t22, -t14, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:12:10
	% EndTime: 2019-12-31 21:12:10
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (138->23), mult. (151->32), div. (0->0), fcn. (174->8), ass. (0->24)
	t22 = sin(qJ(3));
	t24 = cos(qJ(3));
	t38 = pkin(3) + pkin(4);
	t39 = -t22 * qJ(4) - t38 * t24;
	t41 = -pkin(2) + t39;
	t40 = pkin(7) - pkin(8) - r_i_i_C(3);
	t21 = sin(qJ(5));
	t23 = cos(qJ(5));
	t29 = t21 * t24 - t22 * t23;
	t20 = qJ(1) + qJ(2);
	t18 = sin(t20);
	t5 = t29 * t18;
	t28 = t22 * t21 + t23 * t24;
	t6 = t28 * t18;
	t32 = -t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
	t19 = cos(t20);
	t7 = t29 * t19;
	t8 = t28 * t19;
	t31 = -t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
	t30 = -r_i_i_C(1) * t28 + r_i_i_C(2) * t29;
	t27 = qJ(4) * t24 - t38 * t22;
	t26 = t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t40 * t18 - t41 * t19;
	t25 = -t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t41 * t18 + t40 * t19;
	t1 = [-sin(qJ(1)) * pkin(1) + t25, t25, t27 * t19 - t31, t19 * t22, t31; cos(qJ(1)) * pkin(1) + t26, t26, t27 * t18 - t32, t18 * t22, t32; 0, 0, -t30 - t39, -t24, t30;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end