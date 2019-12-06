% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_jacobia_transl_sym_varpar: pkin has to be [2x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:03:48
	% EndTime: 2019-12-05 17:03:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:03:49
	% EndTime: 2019-12-05 17:03:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:03:48
	% EndTime: 2019-12-05 17:03:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(2));
	t1 = sin(qJ(2));
	t3 = [0, -r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; 0, 0, 0, 0, 0; 1, r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:03:48
	% EndTime: 2019-12-05 17:03:49
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (7->5), mult. (20->10), div. (0->0), fcn. (20->4), ass. (0->7)
	t1 = sin(qJ(3));
	t3 = cos(qJ(3));
	t6 = t3 * r_i_i_C(1) - t1 * r_i_i_C(2);
	t5 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t4 = cos(qJ(2));
	t2 = sin(qJ(2));
	t7 = [0, t4 * r_i_i_C(3) - t6 * t2, t5 * t4, 0, 0; 0, 0, -t6, 0, 0; 1, t2 * r_i_i_C(3) + t6 * t4, t5 * t2, 0, 0;];
	Ja_transl = t7;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:03:48
	% EndTime: 2019-12-05 17:03:49
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (31->8), mult. (39->14), div. (0->0), fcn. (39->6), ass. (0->10)
	t4 = qJ(3) + qJ(4);
	t2 = sin(t4);
	t3 = cos(t4);
	t14 = -t3 * r_i_i_C(1) + t2 * r_i_i_C(2);
	t9 = cos(qJ(3)) * pkin(2) - t14;
	t11 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t10 = -sin(qJ(3)) * pkin(2) + t11;
	t8 = cos(qJ(2));
	t6 = sin(qJ(2));
	t1 = [0, t8 * r_i_i_C(3) - t9 * t6, t10 * t8, t11 * t8, 0; 0, 0, -t9, t14, 0; 1, t6 * r_i_i_C(3) + t9 * t8, t10 * t6, t11 * t6, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:03:48
	% EndTime: 2019-12-05 17:03:49
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (62->22), mult. (91->36), div. (0->0), fcn. (99->8), ass. (0->27)
	t12 = qJ(3) + qJ(4);
	t10 = sin(t12);
	t11 = cos(t12);
	t13 = sin(qJ(5));
	t31 = r_i_i_C(2) * t13;
	t35 = r_i_i_C(3) * t11 + t10 * t31;
	t15 = sin(qJ(2));
	t34 = t35 * t15;
	t18 = cos(qJ(2));
	t33 = t35 * t18;
	t16 = cos(qJ(5));
	t32 = r_i_i_C(1) * t16;
	t29 = t10 * r_i_i_C(3);
	t28 = cos(qJ(3)) * pkin(2);
	t27 = t15 * t13;
	t26 = t15 * t16;
	t25 = t18 * t13;
	t24 = t18 * t16;
	t23 = t10 * t32;
	t21 = t28 + t29;
	t20 = -sin(qJ(3)) * pkin(2) - t23;
	t19 = -t29 + (t31 - t32) * t11;
	t4 = t11 * t24 + t27;
	t3 = -t11 * t25 + t26;
	t2 = -t11 * t26 + t25;
	t1 = t11 * t27 + t24;
	t5 = [0, t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t21 * t15, t20 * t18 + t33, -t18 * t23 + t33, t3 * r_i_i_C(1) - t4 * r_i_i_C(2); 0, 0, t19 - t28, t19, (r_i_i_C(1) * t13 + r_i_i_C(2) * t16) * t10; 1, t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t21 * t18, t20 * t15 + t34, -t15 * t23 + t34, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2);];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end