% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:30
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRPR3_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPR3_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR3_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_jacobia_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:38
	% EndTime: 2019-10-24 10:30:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:38
	% EndTime: 2019-10-24 10:30:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:38
	% EndTime: 2019-10-24 10:30:38
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (6->3), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->4)
	t3 = pkin(8) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, -t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0, 0, 0; 0, t2 * r_i_i_C(1) - t1 * r_i_i_C(2), 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:38
	% EndTime: 2019-10-24 10:30:38
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (25->7), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->10)
	t9 = pkin(6) + r_i_i_C(3);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = r_i_i_C(1) * t5 - r_i_i_C(2) * t4;
	t7 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t6 = pkin(2) + t8;
	t3 = pkin(8) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [0, -t1 * t6 + t2 * t9, t7 * t2, 0, 0; 0, t1 * t9 + t2 * t6, t7 * t1, 0, 0; 1, 0, t8, 0, 0;];
	Ja_transl = t10;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:38
	% EndTime: 2019-10-24 10:30:38
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (47->11), mult. (31->12), div. (0->0), fcn. (33->6), ass. (0->11)
	t8 = qJ(3) + pkin(9);
	t3 = sin(t8);
	t5 = cos(t8);
	t15 = t5 * r_i_i_C(1) - t3 * r_i_i_C(2) + cos(qJ(3)) * pkin(3);
	t14 = r_i_i_C(3) + qJ(4) + pkin(6);
	t12 = pkin(2) + t15;
	t11 = -sin(qJ(3)) * pkin(3) - r_i_i_C(1) * t3 - r_i_i_C(2) * t5;
	t7 = pkin(8) + qJ(2);
	t4 = cos(t7);
	t2 = sin(t7);
	t1 = [0, -t12 * t2 + t14 * t4, t11 * t4, t2, 0; 0, t12 * t4 + t14 * t2, t11 * t2, -t4, 0; 1, 0, t15, 0, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:30:38
	% EndTime: 2019-10-24 10:30:38
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (88->14), mult. (46->16), div. (0->0), fcn. (48->8), ass. (0->14)
	t13 = qJ(3) + pkin(9);
	t9 = qJ(5) + t13;
	t5 = sin(t9);
	t6 = cos(t9);
	t17 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
	t21 = t17 + cos(qJ(3)) * pkin(3) + pkin(4) * cos(t13);
	t19 = r_i_i_C(3) + pkin(7) + qJ(4) + pkin(6);
	t16 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
	t15 = pkin(2) + t21;
	t14 = t16 - pkin(4) * sin(t13) - sin(qJ(3)) * pkin(3);
	t12 = pkin(8) + qJ(2);
	t8 = cos(t12);
	t7 = sin(t12);
	t1 = [0, -t15 * t7 + t19 * t8, t14 * t8, t7, t16 * t8; 0, t15 * t8 + t19 * t7, t14 * t7, -t8, t16 * t7; 1, 0, t21, 0, t17;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end