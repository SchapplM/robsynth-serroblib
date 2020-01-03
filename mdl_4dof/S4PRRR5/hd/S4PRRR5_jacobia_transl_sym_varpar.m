% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PRRR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4PRRR5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRR5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->2), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t3 = sin(qJ(2));
	t4 = cos(qJ(2));
	t5 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t1 = [0, t5 * cos(pkin(7)), 0, 0; 0, t5 * sin(pkin(7)), 0, 0; 1, t4 * r_i_i_C(1) - t3 * r_i_i_C(2), 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (21->5), mult. (25->10), div. (0->0), fcn. (25->6), ass. (0->9)
	t4 = qJ(2) + qJ(3);
	t2 = sin(t4);
	t3 = cos(t4);
	t10 = r_i_i_C(1) * t3 - r_i_i_C(2) * t2;
	t9 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
	t8 = -sin(qJ(2)) * pkin(2) + t9;
	t6 = cos(pkin(7));
	t5 = sin(pkin(7));
	t1 = [0, t8 * t6, t9 * t6, 0; 0, t8 * t5, t9 * t5, 0; 1, cos(qJ(2)) * pkin(2) + t10, t10, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:33:55
	% EndTime: 2019-12-31 16:33:55
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (70->19), mult. (87->30), div. (0->0), fcn. (91->8), ass. (0->21)
	t13 = qJ(2) + qJ(3);
	t11 = sin(t13);
	t12 = cos(t13);
	t16 = sin(qJ(4));
	t32 = r_i_i_C(2) * t16;
	t35 = pkin(6) + r_i_i_C(3);
	t36 = t11 * t32 + t12 * t35;
	t18 = cos(qJ(4));
	t34 = -t18 * r_i_i_C(1) - pkin(3);
	t14 = sin(pkin(7));
	t28 = t14 * t16;
	t27 = t14 * t18;
	t15 = cos(pkin(7));
	t26 = t15 * t16;
	t25 = t15 * t18;
	t24 = t36 * t14;
	t23 = t36 * t15;
	t21 = t34 * t11;
	t20 = t35 * t11 + (-t32 - t34) * t12;
	t19 = -sin(qJ(2)) * pkin(2) + t21;
	t1 = [0, t19 * t15 + t23, t15 * t21 + t23, (-t12 * t26 + t27) * r_i_i_C(1) + (-t12 * t25 - t28) * r_i_i_C(2); 0, t19 * t14 + t24, t14 * t21 + t24, (-t12 * t28 - t25) * r_i_i_C(1) + (-t12 * t27 + t26) * r_i_i_C(2); 1, cos(qJ(2)) * pkin(2) + t20, t20, (-r_i_i_C(1) * t16 - r_i_i_C(2) * t18) * t11;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,4);
end