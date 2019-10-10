% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RPPP1
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
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RPPP1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPP1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPP1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobia_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (13->12), mult. (28->19), div. (0->0), fcn. (35->6), ass. (0->12)
	t1 = sin(pkin(6));
	t5 = sin(qJ(1));
	t11 = t5 * t1;
	t3 = cos(pkin(6));
	t10 = t5 * t3;
	t6 = cos(qJ(1));
	t9 = t6 * t1;
	t8 = t6 * t3;
	t2 = sin(pkin(4));
	t7 = t2 * (r_i_i_C(3) + qJ(2));
	t4 = cos(pkin(4));
	t12 = [(-t4 * t9 - t10) * r_i_i_C(1) + (-t4 * t8 + t11) * r_i_i_C(2) - t5 * pkin(1) + t6 * t7, t5 * t2, 0, 0; (-t11 * t4 + t8) * r_i_i_C(1) + (-t10 * t4 - t9) * r_i_i_C(2) + t6 * pkin(1) + t5 * t7, -t6 * t2, 0, 0; 0, t4, 0, 0;];
	Ja_transl = t12;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (24->15), mult. (51->21), div. (0->0), fcn. (67->6), ass. (0->15)
	t5 = sin(pkin(6));
	t9 = sin(qJ(1));
	t16 = t9 * t5;
	t7 = cos(pkin(6));
	t15 = t9 * t7;
	t14 = pkin(2) - r_i_i_C(2);
	t10 = cos(qJ(1));
	t8 = cos(pkin(4));
	t13 = t10 * t8;
	t12 = r_i_i_C(3) + qJ(3);
	t6 = sin(pkin(4));
	t11 = (r_i_i_C(1) + qJ(2)) * t6;
	t3 = t10 * t5 + t15 * t8;
	t1 = -t13 * t7 + t16;
	t2 = [-t9 * pkin(1) - t14 * (t5 * t13 + t15) + t10 * t11 - t12 * t1, t9 * t6, t3, 0; t10 * pkin(1) + t9 * t11 + t14 * (t10 * t7 - t8 * t16) + t12 * t3, -t10 * t6, t1, 0; 0, t8, -t6 * t7, 0;];
	Ja_transl = t2;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (32->16), mult. (70->22), div. (0->0), fcn. (93->6), ass. (0->17)
	t5 = sin(pkin(6));
	t9 = sin(qJ(1));
	t16 = t9 * t5;
	t7 = cos(pkin(6));
	t15 = t9 * t7;
	t10 = cos(qJ(1));
	t8 = cos(pkin(4));
	t14 = t10 * t8;
	t13 = r_i_i_C(2) + qJ(3);
	t12 = pkin(2) + r_i_i_C(3) + qJ(4);
	t6 = sin(pkin(4));
	t11 = (pkin(3) + r_i_i_C(1) + qJ(2)) * t6;
	t4 = t10 * t7 - t8 * t16;
	t3 = t10 * t5 + t8 * t15;
	t2 = t5 * t14 + t15;
	t1 = -t7 * t14 + t16;
	t17 = [-t9 * pkin(1) - t13 * t1 + t10 * t11 - t12 * t2, t9 * t6, t3, t4; t10 * pkin(1) + t9 * t11 + t12 * t4 + t13 * t3, -t10 * t6, t1, t2; 0, t8, -t6 * t7, t6 * t5;];
	Ja_transl = t17;
else
	Ja_transl=NaN(3,4);
end