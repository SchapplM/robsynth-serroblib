% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RRPR5
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RRPR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRPR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_jacobia_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:03:41
	% EndTime: 2019-12-31 17:03:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:03:41
	% EndTime: 2019-12-31 17:03:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:03:41
	% EndTime: 2019-12-31 17:03:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t4 = qJ(1) + qJ(2);
	t2 = sin(t4);
	t3 = cos(t4);
	t6 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t5 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t1 = [-sin(qJ(1)) * pkin(1) + t5, t5, 0, 0; cos(qJ(1)) * pkin(1) + t6, t6, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:03:41
	% EndTime: 2019-12-31 17:03:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (33->8), mult. (18->6), div. (0->0), fcn. (20->4), ass. (0->8)
	t13 = pkin(2) - r_i_i_C(2);
	t12 = r_i_i_C(3) + qJ(3);
	t9 = qJ(1) + qJ(2);
	t7 = sin(t9);
	t8 = cos(t9);
	t11 = t12 * t8 - t13 * t7;
	t10 = t12 * t7 + t13 * t8;
	t1 = [-sin(qJ(1)) * pkin(1) + t11, t11, t7, 0; cos(qJ(1)) * pkin(1) + t10, t10, -t8, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 17:03:41
	% EndTime: 2019-12-31 17:03:41
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (56->11), mult. (44->12), div. (0->0), fcn. (46->6), ass. (0->12)
	t13 = sin(qJ(4));
	t14 = cos(qJ(4));
	t22 = -t13 * r_i_i_C(1) - t14 * r_i_i_C(2);
	t21 = pkin(2) + pkin(6) + r_i_i_C(3);
	t20 = qJ(3) - t22;
	t12 = qJ(1) + qJ(2);
	t10 = sin(t12);
	t11 = cos(t12);
	t17 = t20 * t10 + t21 * t11;
	t16 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13;
	t15 = -t21 * t10 + t20 * t11;
	t1 = [-sin(qJ(1)) * pkin(1) + t15, t15, t10, t16 * t10; cos(qJ(1)) * pkin(1) + t17, t17, -t11, -t16 * t11; 0, 0, 0, t22;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,4);
end