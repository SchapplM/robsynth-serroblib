% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RPPR2
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RPPR2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPR2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPR2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_jacobia_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:41:54
	% EndTime: 2019-10-09 20:41:54
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:41:54
	% EndTime: 2019-10-09 20:41:54
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
	% StartTime: 2019-10-09 20:41:54
	% EndTime: 2019-10-09 20:41:54
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (7->5), mult. (8->4), div. (0->0), fcn. (10->2), ass. (0->5)
	t4 = pkin(1) + r_i_i_C(1);
	t3 = r_i_i_C(3) + qJ(2);
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t5 = [-t1 * t4 + t2 * t3, t1, 0, 0; t1 * t3 + t2 * t4, -t2, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t5;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:41:54
	% EndTime: 2019-10-09 20:41:54
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (13->10), mult. (18->12), div. (0->0), fcn. (24->4), ass. (0->8)
	t7 = pkin(1) + pkin(2);
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t4 = cos(pkin(6));
	t3 = sin(pkin(6));
	t2 = t3 * t5 + t4 * t6;
	t1 = t3 * t6 - t4 * t5;
	t8 = [t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + t6 * qJ(2) - t5 * t7, t5, 0, 0; t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + t5 * qJ(2) + t6 * t7, -t6, 0, 0; 0, 0, -1, 0;];
	Ja_transl = t8;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:41:54
	% EndTime: 2019-10-09 20:41:54
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (39->13), mult. (36->14), div. (0->0), fcn. (46->6), ass. (0->12)
	t16 = pkin(1) + cos(pkin(6)) * pkin(3) + pkin(2);
	t15 = pkin(6) + qJ(4);
	t14 = pkin(3) * sin(pkin(6)) + qJ(2);
	t13 = cos(t15);
	t12 = sin(t15);
	t8 = sin(qJ(1));
	t9 = cos(qJ(1));
	t1 = -t8 * t12 - t9 * t13;
	t2 = t9 * t12 - t8 * t13;
	t11 = -t2 * r_i_i_C(1) + t1 * r_i_i_C(2);
	t10 = t1 * r_i_i_C(1) + t2 * r_i_i_C(2);
	t3 = [t14 * t9 - t16 * t8 - t11, t8, 0, t11; t14 * t8 + t16 * t9 - t10, -t9, 0, t10; 0, 0, -1, 0;];
	Ja_transl = t3;
else
	Ja_transl=NaN(3,4);
end