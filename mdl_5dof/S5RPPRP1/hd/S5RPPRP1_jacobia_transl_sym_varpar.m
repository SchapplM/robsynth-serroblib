% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPPRP1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRP1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRP1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_jacobia_transl_sym_varpar: pkin has to be [8x1] (double)');
Ja_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
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
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->5), mult. (6->6), div. (0->0), fcn. (6->4), ass. (0->4)
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-t1 * r_i_i_C(1) - t2 * r_i_i_C(2) - sin(qJ(1)) * pkin(1), 0, 0, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + cos(qJ(1)) * pkin(1), 0, 0, 0, 0; 0, 1, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (23->9), mult. (16->8), div. (0->0), fcn. (18->6), ass. (0->6)
	t7 = r_i_i_C(3) + qJ(3);
	t6 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(2);
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [-sin(qJ(1)) * pkin(1) + t7 * t2 - t6 * t1, 0, t1, 0, 0; cos(qJ(1)) * pkin(1) + t7 * t1 + t6 * t2, 0, -t2, 0, 0; 0, 1, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (54->20), mult. (54->27), div. (0->0), fcn. (64->8), ass. (0->17)
	t8 = sin(pkin(8));
	t9 = cos(pkin(8));
	t18 = pkin(3) * t9 + pkin(2) + (pkin(6) + r_i_i_C(3)) * t8;
	t10 = sin(qJ(4));
	t7 = qJ(1) + pkin(7);
	t5 = sin(t7);
	t16 = t5 * t10;
	t11 = cos(qJ(4));
	t15 = t5 * t11;
	t6 = cos(t7);
	t14 = t6 * t10;
	t13 = t6 * t11;
	t4 = t9 * t13 + t16;
	t3 = -t9 * t14 + t15;
	t2 = -t9 * t15 + t14;
	t1 = t9 * t16 + t13;
	t12 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t6 * qJ(3) - t18 * t5, 0, t5, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t5 * qJ(3) + t18 * t6, 0, -t6, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0; 0, 1, 0, (-r_i_i_C(1) * t10 - r_i_i_C(2) * t11) * t8, 0;];
	Ja_transl = t12;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:13:04
	% EndTime: 2022-01-23 09:13:04
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (74->25), mult. (74->33), div. (0->0), fcn. (87->8), ass. (0->17)
	t20 = pkin(4) + r_i_i_C(1);
	t10 = cos(pkin(8));
	t13 = cos(qJ(4));
	t9 = sin(pkin(8));
	t19 = (r_i_i_C(3) + qJ(5) + pkin(6)) * t9 + t10 * (pkin(4) * t13 + pkin(3)) + pkin(2);
	t12 = sin(qJ(4));
	t17 = t10 * t12;
	t16 = t10 * t13;
	t14 = pkin(4) * t12 + qJ(3);
	t8 = qJ(1) + pkin(7);
	t6 = sin(t8);
	t7 = cos(t8);
	t3 = t13 * t6 - t7 * t17;
	t1 = t13 * t7 + t6 * t17;
	t4 = t6 * t12 + t7 * t16;
	t2 = t7 * t12 - t6 * t16;
	t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t14 * t7 - t19 * t6, 0, t6, -t4 * r_i_i_C(2) + t20 * t3, t7 * t9; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t14 * t6 + t19 * t7, 0, -t7, t2 * r_i_i_C(2) - t20 * t1, t6 * t9; 0, 1, 0, (-r_i_i_C(2) * t13 - t20 * t12) * t9, -t10;];
	Ja_transl = t5;
end