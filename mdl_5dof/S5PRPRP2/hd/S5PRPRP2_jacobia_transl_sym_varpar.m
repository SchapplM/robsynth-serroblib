% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRPRP2
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRPRP2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRPRP2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:31:30
	% EndTime: 2019-12-05 15:31:30
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:31:30
	% EndTime: 2019-12-05 15:31:30
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:31:30
	% EndTime: 2019-12-05 15:31:30
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (6->3), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->4)
	t3 = pkin(7) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, -t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0, 0; 0, t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:31:30
	% EndTime: 2019-12-05 15:31:30
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (21->7), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->6)
	t7 = r_i_i_C(3) + qJ(3);
	t6 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(2);
	t3 = pkin(7) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, -t6 * t1 + t7 * t2, t1, 0, 0; 0, t7 * t1 + t6 * t2, -t2, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:31:30
	% EndTime: 2019-12-05 15:31:30
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (52->18), mult. (52->25), div. (0->0), fcn. (62->6), ass. (0->17)
	t8 = sin(pkin(8));
	t9 = cos(pkin(8));
	t18 = pkin(3) * t9 + pkin(2) + (pkin(6) + r_i_i_C(3)) * t8;
	t10 = sin(qJ(4));
	t7 = pkin(7) + qJ(2);
	t5 = sin(t7);
	t16 = t5 * t10;
	t11 = cos(qJ(4));
	t15 = t5 * t11;
	t6 = cos(t7);
	t14 = t6 * t10;
	t13 = t6 * t11;
	t4 = t13 * t9 + t16;
	t3 = -t14 * t9 + t15;
	t2 = -t15 * t9 + t14;
	t1 = t16 * t9 + t13;
	t12 = [0, t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t6 * qJ(3) - t18 * t5, t5, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0; 0, t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t5 * qJ(3) + t18 * t6, -t6, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0; 1, 0, 0, (-r_i_i_C(1) * t10 - r_i_i_C(2) * t11) * t8, 0;];
	Ja_transl = t12;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:31:30
	% EndTime: 2019-12-05 15:31:30
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (72->23), mult. (72->31), div. (0->0), fcn. (85->6), ass. (0->17)
	t20 = pkin(4) + r_i_i_C(1);
	t10 = cos(pkin(8));
	t13 = cos(qJ(4));
	t9 = sin(pkin(8));
	t19 = (r_i_i_C(3) + qJ(5) + pkin(6)) * t9 + t10 * (pkin(4) * t13 + pkin(3)) + pkin(2);
	t12 = sin(qJ(4));
	t17 = t10 * t12;
	t16 = t10 * t13;
	t14 = pkin(4) * t12 + qJ(3);
	t8 = pkin(7) + qJ(2);
	t6 = sin(t8);
	t7 = cos(t8);
	t3 = t13 * t6 - t7 * t17;
	t1 = t13 * t7 + t6 * t17;
	t4 = t6 * t12 + t7 * t16;
	t2 = t7 * t12 - t6 * t16;
	t5 = [0, t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t14 * t7 - t19 * t6, t6, -t4 * r_i_i_C(2) + t20 * t3, t7 * t9; 0, t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t14 * t6 + t19 * t7, -t7, t2 * r_i_i_C(2) - t1 * t20, t6 * t9; 1, 0, 0, (-r_i_i_C(2) * t13 - t12 * t20) * t9, -t10;];
	Ja_transl = t5;
else
	Ja_transl=NaN(3,5);
end