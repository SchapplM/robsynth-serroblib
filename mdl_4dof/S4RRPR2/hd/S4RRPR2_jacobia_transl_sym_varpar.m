% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RRPR2
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:51
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RRPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_jacobia_transl_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:50
	% EndTime: 2019-10-09 20:51:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:50
	% EndTime: 2019-10-09 20:51:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:50
	% EndTime: 2019-10-09 20:51:50
	% DurationCPUTime: 0.08s
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
	% StartTime: 2019-10-09 20:51:50
	% EndTime: 2019-10-09 20:51:50
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (33->8), mult. (18->6), div. (0->0), fcn. (20->4), ass. (0->8)
	t13 = pkin(2) + r_i_i_C(1);
	t12 = qJ(3) + r_i_i_C(3);
	t9 = qJ(1) + qJ(2);
	t7 = sin(t9);
	t8 = cos(t9);
	t11 = t12 * t7 + t13 * t8;
	t10 = t12 * t8 - t13 * t7;
	t1 = [-sin(qJ(1)) * pkin(1) + t10, t10, t7, 0; cos(qJ(1)) * pkin(1) + t11, t11, -t8, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:50
	% EndTime: 2019-10-09 20:51:50
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (71->13), mult. (50->14), div. (0->0), fcn. (64->6), ass. (0->13)
	t23 = pkin(2) + pkin(3);
	t16 = qJ(1) + qJ(2);
	t14 = sin(t16);
	t15 = cos(t16);
	t19 = sin(qJ(4));
	t20 = cos(qJ(4));
	t5 = -t14 * t19 - t15 * t20;
	t6 = -t14 * t20 + t15 * t19;
	t22 = -t6 * r_i_i_C(1) + t5 * r_i_i_C(2);
	t21 = t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
	t18 = t14 * qJ(3) + t23 * t15 - t21;
	t17 = t15 * qJ(3) - t23 * t14 - t22;
	t1 = [-sin(qJ(1)) * pkin(1) + t17, t17, t14, t22; cos(qJ(1)) * pkin(1) + t18, t18, -t15, t21; 0, 0, 0, 0;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,4);
end