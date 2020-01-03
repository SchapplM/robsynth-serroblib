% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRPPR4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPPR4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPR4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:28:04
	% EndTime: 2019-12-31 19:28:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:28:04
	% EndTime: 2019-12-31 19:28:04
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
	% StartTime: 2019-12-31 19:28:04
	% EndTime: 2019-12-31 19:28:04
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
	% StartTime: 2019-12-31 19:28:04
	% EndTime: 2019-12-31 19:28:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (33->8), mult. (18->6), div. (0->0), fcn. (20->4), ass. (0->8)
	t13 = pkin(2) + r_i_i_C(1);
	t12 = qJ(3) + r_i_i_C(3);
	t9 = qJ(1) + qJ(2);
	t7 = sin(t9);
	t8 = cos(t9);
	t11 = t12 * t7 + t13 * t8;
	t10 = t12 * t8 - t13 * t7;
	t1 = [-sin(qJ(1)) * pkin(1) + t10, t10, t7, 0, 0; cos(qJ(1)) * pkin(1) + t11, t11, -t8, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:28:04
	% EndTime: 2019-12-31 19:28:04
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (57->13), mult. (38->14), div. (0->0), fcn. (48->6), ass. (0->11)
	t18 = pkin(2) + pkin(3);
	t13 = qJ(1) + qJ(2);
	t11 = sin(t13);
	t12 = cos(t13);
	t14 = sin(pkin(8));
	t15 = cos(pkin(8));
	t5 = -t11 * t15 + t12 * t14;
	t6 = t11 * t14 + t12 * t15;
	t17 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2) + t11 * qJ(3) + t18 * t12;
	t16 = t5 * r_i_i_C(1) + t6 * r_i_i_C(2) + t12 * qJ(3) - t18 * t11;
	t1 = [-sin(qJ(1)) * pkin(1) + t16, t16, t11, 0, 0; cos(qJ(1)) * pkin(1) + t17, t17, -t12, 0, 0; 0, 0, 0, -1, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:28:04
	% EndTime: 2019-12-31 19:28:04
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (120->17), mult. (100->20), div. (0->0), fcn. (126->8), ass. (0->17)
	t31 = pkin(2) + pkin(3);
	t30 = pkin(7) + r_i_i_C(3);
	t16 = sin(qJ(5));
	t17 = cos(qJ(5));
	t29 = -t17 * r_i_i_C(1) + t16 * r_i_i_C(2);
	t28 = pkin(4) - t29;
	t25 = cos(pkin(8));
	t24 = sin(pkin(8));
	t23 = qJ(1) + qJ(2);
	t22 = cos(t23);
	t21 = sin(t23);
	t20 = r_i_i_C(1) * t16 + r_i_i_C(2) * t17;
	t7 = -t21 * t24 - t22 * t25;
	t8 = -t21 * t25 + t22 * t24;
	t19 = t21 * qJ(3) + t31 * t22 - t28 * t7 + t30 * t8;
	t18 = t22 * qJ(3) - t31 * t21 + t28 * t8 + t30 * t7;
	t1 = [-sin(qJ(1)) * pkin(1) + t18, t18, t21, 0, t20 * t7; cos(qJ(1)) * pkin(1) + t19, t19, -t22, 0, t20 * t8; 0, 0, 0, -1, t29;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end