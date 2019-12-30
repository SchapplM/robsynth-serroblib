% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:37
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5PRRPP4_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPP4_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP4_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:37:36
	% EndTime: 2019-12-29 15:37:36
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:37:36
	% EndTime: 2019-12-29 15:37:36
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:37:36
	% EndTime: 2019-12-29 15:37:36
	% DurationCPUTime: 0.11s
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
	% StartTime: 2019-12-29 15:37:36
	% EndTime: 2019-12-29 15:37:36
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (25->7), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->10)
	t9 = pkin(6) + r_i_i_C(3);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t7 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t6 = pkin(2) + t8;
	t3 = pkin(7) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [0, -t6 * t1 + t9 * t2, t7 * t2, 0, 0; 0, t9 * t1 + t6 * t2, t7 * t1, 0, 0; 1, 0, t8, 0, 0;];
	Ja_transl = t10;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:37:31
	% EndTime: 2019-12-29 15:37:31
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (46->10), mult. (44->12), div. (0->0), fcn. (47->4), ass. (0->12)
	t10 = pkin(3) + r_i_i_C(1);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = r_i_i_C(3) + qJ(4);
	t7 = t10 * t5 + t8 * t4;
	t11 = pkin(2) + t7;
	t9 = pkin(6) + r_i_i_C(2);
	t6 = -t10 * t4 + t8 * t5;
	t3 = pkin(7) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t12 = [0, -t11 * t1 + t9 * t2, t6 * t2, t2 * t4, 0; 0, t9 * t1 + t11 * t2, t6 * t1, t1 * t4, 0; 1, 0, t7, -t5, 0;];
	Ja_transl = t12;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:37:36
	% EndTime: 2019-12-29 15:37:36
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (62->11), mult. (55->12), div. (0->0), fcn. (60->4), ass. (0->12)
	t10 = r_i_i_C(2) + qJ(4);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t9 = pkin(3) + pkin(4) + r_i_i_C(1);
	t7 = t10 * t4 + t9 * t5;
	t11 = pkin(2) + t7;
	t8 = pkin(6) - r_i_i_C(3) - qJ(5);
	t6 = t10 * t5 - t9 * t4;
	t3 = pkin(7) + qJ(2);
	t2 = cos(t3);
	t1 = sin(t3);
	t12 = [0, -t11 * t1 + t8 * t2, t6 * t2, t2 * t4, -t1; 0, t8 * t1 + t11 * t2, t6 * t1, t1 * t4, t2; 1, 0, t7, -t5, 0;];
	Ja_transl = t12;
else
	Ja_transl=NaN(3,5);
end