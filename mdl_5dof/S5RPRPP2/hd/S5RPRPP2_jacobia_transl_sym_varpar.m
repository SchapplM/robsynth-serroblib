% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 16:35
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRPP2_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPP2_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPP2_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:35:18
	% EndTime: 2019-12-29 16:35:18
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:35:23
	% EndTime: 2019-12-29 16:35:23
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:35:23
	% EndTime: 2019-12-29 16:35:23
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-12-29 16:35:18
	% EndTime: 2019-12-29 16:35:18
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (27->9), mult. (26->12), div. (0->0), fcn. (26->6), ass. (0->10)
	t9 = pkin(6) + r_i_i_C(3);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
	t7 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
	t6 = pkin(2) + t8;
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t10 = [-sin(qJ(1)) * pkin(1) + t9 * t2 - t6 * t1, 0, t7 * t2, 0, 0; cos(qJ(1)) * pkin(1) + t9 * t1 + t6 * t2, 0, t7 * t1, 0, 0; 0, 1, t8, 0, 0;];
	Ja_transl = t10;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:35:23
	% EndTime: 2019-12-29 16:35:23
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (48->12), mult. (46->14), div. (0->0), fcn. (49->6), ass. (0->12)
	t10 = pkin(3) + r_i_i_C(1);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t8 = r_i_i_C(3) + qJ(4);
	t7 = t10 * t5 + t8 * t4;
	t11 = pkin(2) + t7;
	t9 = pkin(6) + r_i_i_C(2);
	t6 = -t10 * t4 + t8 * t5;
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t12 = [-sin(qJ(1)) * pkin(1) + t9 * t2 - t11 * t1, 0, t6 * t2, t2 * t4, 0; cos(qJ(1)) * pkin(1) + t9 * t1 + t11 * t2, 0, t6 * t1, t1 * t4, 0; 0, 1, t7, -t5, 0;];
	Ja_transl = t12;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 16:35:23
	% EndTime: 2019-12-29 16:35:23
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (64->13), mult. (57->14), div. (0->0), fcn. (62->6), ass. (0->12)
	t10 = r_i_i_C(2) + qJ(4);
	t4 = sin(qJ(3));
	t5 = cos(qJ(3));
	t9 = pkin(3) + pkin(4) + r_i_i_C(1);
	t7 = t10 * t4 + t9 * t5;
	t11 = pkin(2) + t7;
	t8 = pkin(6) - r_i_i_C(3) - qJ(5);
	t6 = t10 * t5 - t9 * t4;
	t3 = qJ(1) + pkin(7);
	t2 = cos(t3);
	t1 = sin(t3);
	t12 = [-sin(qJ(1)) * pkin(1) + t8 * t2 - t11 * t1, 0, t6 * t2, t2 * t4, -t1; cos(qJ(1)) * pkin(1) + t8 * t1 + t11 * t2, 0, t6 * t1, t1 * t4, t2; 0, 1, t7, -t5, 0;];
	Ja_transl = t12;
else
	Ja_transl=NaN(3,5);
end