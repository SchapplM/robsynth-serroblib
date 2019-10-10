% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S3RRR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% 
% Output:
% Ja_transl [3x3]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S3RRR1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_jacobia_transl_sym_varpar: qJ has to be [3x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S3RRR1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S3RRR1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_jacobia_transl_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:13:34
	% EndTime: 2019-10-09 20:13:34
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0; 0, 0, 0; 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:13:34
	% EndTime: 2019-10-09 20:13:34
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0; 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:13:34
	% EndTime: 2019-10-09 20:13:34
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (14->5), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->6)
	t4 = qJ(1) + qJ(2);
	t2 = sin(t4);
	t3 = cos(t4);
	t6 = t3 * r_i_i_C(1) - t2 * r_i_i_C(2);
	t5 = -t2 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t1 = [-sin(qJ(1)) * pkin(1) + t5, t5, 0; cos(qJ(1)) * pkin(1) + t6, t6, 0; 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:13:34
	% EndTime: 2019-10-09 20:13:34
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (40->8), mult. (18->8), div. (0->0), fcn. (18->6), ass. (0->9)
	t7 = qJ(1) + qJ(2);
	t6 = qJ(3) + t7;
	t3 = sin(t6);
	t4 = cos(t6);
	t11 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t10 = t11 + pkin(2) * cos(t7);
	t9 = -t3 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t8 = -pkin(2) * sin(t7) + t9;
	t1 = [-sin(qJ(1)) * pkin(1) + t8, t8, t9; cos(qJ(1)) * pkin(1) + t10, t10, t11; 0, 0, 0;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,3);
end