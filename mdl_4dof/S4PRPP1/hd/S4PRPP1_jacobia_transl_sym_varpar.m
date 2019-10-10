% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4PRPP1
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
%   pkin=[a2,a3,a4,d2,theta1]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:26
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4PRPP1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRPP1_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRPP1_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_jacobia_transl_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:26:20
	% EndTime: 2019-10-09 20:26:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:26:20
	% EndTime: 2019-10-09 20:26:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:26:20
	% EndTime: 2019-10-09 20:26:20
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (6->3), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->4)
	t3 = pkin(5) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t4 = [0, -r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0; 0, r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0; 1, 0, 0, 0;];
	Ja_transl = t4;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:26:20
	% EndTime: 2019-10-09 20:26:20
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (17->6), mult. (8->4), div. (0->0), fcn. (10->2), ass. (0->6)
	t5 = pkin(2) - r_i_i_C(2);
	t4 = r_i_i_C(3) + qJ(3);
	t3 = pkin(5) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t6 = [0, -t5 * t1 + t4 * t2, t1, 0; 0, t4 * t1 + t5 * t2, -t2, 0; 1, 0, 0, 0;];
	Ja_transl = t6;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:26:20
	% EndTime: 2019-10-09 20:26:20
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (23->6), mult. (10->4), div. (0->0), fcn. (14->2), ass. (0->6)
	t5 = r_i_i_C(2) + qJ(3);
	t4 = pkin(2) + r_i_i_C(3) + qJ(4);
	t3 = pkin(5) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t6 = [0, -t4 * t1 + t5 * t2, t1, t2; 0, t5 * t1 + t4 * t2, -t2, t1; 1, 0, 0, 0;];
	Ja_transl = t6;
else
	Ja_transl=NaN(3,4);
end