% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR12
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRPR12_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR12_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR12_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:31:04
	% EndTime: 2019-12-31 18:31:04
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:31:04
	% EndTime: 2019-12-31 18:31:04
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
	% StartTime: 2019-12-31 18:31:04
	% EndTime: 2019-12-31 18:31:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t3 * t5 + t4 * t6, t3, 0, 0, 0; t3 * t6 + t4 * t5, -t4, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:31:04
	% EndTime: 2019-12-31 18:31:04
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (26->9), mult. (26->11), div. (0->0), fcn. (28->5), ass. (0->10)
	t11 = r_i_i_C(3) + pkin(6) + qJ(2);
	t4 = pkin(8) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = cos(pkin(8)) * pkin(2) + pkin(1) + t10;
	t7 = cos(qJ(1));
	t6 = sin(qJ(1));
	t1 = [t11 * t7 - t8 * t6, t6, t9 * t7, 0, 0; t11 * t6 + t8 * t7, -t7, t9 * t6, 0, 0; 0, 0, t10, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:31:04
	% EndTime: 2019-12-31 18:31:04
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (62->15), mult. (67->17), div. (0->0), fcn. (76->7), ass. (0->14)
	t5 = sin(pkin(9));
	t6 = cos(pkin(9));
	t13 = r_i_i_C(1) * t6 - r_i_i_C(2) * t5 + pkin(3);
	t14 = r_i_i_C(3) + qJ(4);
	t4 = pkin(8) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t11 = t13 * t3 + t14 * t2;
	t15 = cos(pkin(8)) * pkin(2) + pkin(1) + t11;
	t12 = r_i_i_C(1) * t5 + r_i_i_C(2) * t6 + pkin(6) + qJ(2);
	t10 = -t13 * t2 + t14 * t3;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [t12 * t9 - t15 * t8, t8, t10 * t9, t9 * t2, 0; t12 * t8 + t15 * t9, -t9, t10 * t8, t8 * t2, 0; 0, 0, t11, -t3, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:31:04
	% EndTime: 2019-12-31 18:31:04
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (111->28), mult. (96->39), div. (0->0), fcn. (109->9), ass. (0->22)
	t12 = pkin(8) + qJ(3);
	t10 = cos(t12);
	t25 = r_i_i_C(3) + pkin(7) + qJ(4);
	t8 = sin(t12);
	t22 = t25 * t8;
	t5 = cos(pkin(9)) * pkin(4) + pkin(3);
	t26 = t22 + t10 * t5 + cos(pkin(8)) * pkin(2) + pkin(1);
	t16 = sin(qJ(1));
	t24 = t10 * t16;
	t17 = cos(qJ(1));
	t23 = t10 * t17;
	t20 = pkin(4) * sin(pkin(9)) + pkin(6) + qJ(2);
	t11 = pkin(9) + qJ(5);
	t7 = sin(t11);
	t9 = cos(t11);
	t19 = r_i_i_C(1) * t9 - r_i_i_C(2) * t7 + t5;
	t18 = t25 * t10 - t19 * t8;
	t4 = t16 * t7 + t9 * t23;
	t3 = t16 * t9 - t7 * t23;
	t2 = t17 * t7 - t9 * t24;
	t1 = t17 * t9 + t7 * t24;
	t6 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t26 * t16 + t20 * t17, t16, t18 * t17, t17 * t8, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t20 * t16 + t26 * t17, -t17, t18 * t16, t16 * t8, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 0, t19 * t10 + t22, -t10, (-r_i_i_C(1) * t7 - r_i_i_C(2) * t9) * t8;];
	Ja_transl = t6;
else
	Ja_transl=NaN(3,5);
end