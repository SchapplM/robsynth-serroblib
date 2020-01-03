% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RPPR6
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% Ja_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RPPR6_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_jacobia_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPR6_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPR6_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_jacobia_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:40:54
	% EndTime: 2019-12-31 16:40:54
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:40:54
	% EndTime: 2019-12-31 16:40:54
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
	% StartTime: 2019-12-31 16:40:54
	% EndTime: 2019-12-31 16:40:54
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->6), mult. (14->6), div. (0->0), fcn. (16->4), ass. (0->5)
	t6 = r_i_i_C(3) + qJ(2);
	t5 = r_i_i_C(1) * cos(pkin(6)) - r_i_i_C(2) * sin(pkin(6)) + pkin(1);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t1 = [-t5 * t3 + t6 * t4, t3, 0, 0; t6 * t3 + t5 * t4, -t4, 0, 0; 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:40:54
	% EndTime: 2019-12-31 16:40:54
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (14->9), mult. (24->8), div. (0->0), fcn. (29->4), ass. (0->7)
	t1 = sin(pkin(6));
	t2 = cos(pkin(6));
	t8 = pkin(1) + (pkin(2) + r_i_i_C(1)) * t2 + (r_i_i_C(3) + qJ(3)) * t1;
	t6 = r_i_i_C(2) + qJ(2);
	t4 = cos(qJ(1));
	t3 = sin(qJ(1));
	t5 = [-t8 * t3 + t6 * t4, t3, t4 * t1, 0; t6 * t3 + t8 * t4, -t4, t3 * t1, 0; 0, 0, -t2, 0;];
	Ja_transl = t5;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 16:40:54
	% EndTime: 2019-12-31 16:40:54
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (31->17), mult. (68->26), div. (0->0), fcn. (83->6), ass. (0->15)
	t5 = sin(pkin(6));
	t6 = cos(pkin(6));
	t16 = qJ(3) * t5 + pkin(1) + (pkin(2) + pkin(3)) * t6;
	t14 = -pkin(5) - r_i_i_C(3) + qJ(2);
	t7 = sin(qJ(4));
	t9 = cos(qJ(4));
	t12 = t5 * t9 - t6 * t7;
	t11 = t5 * t7 + t6 * t9;
	t10 = cos(qJ(1));
	t8 = sin(qJ(1));
	t4 = t11 * t10;
	t3 = t12 * t10;
	t2 = t11 * t8;
	t1 = t12 * t8;
	t13 = [-t2 * r_i_i_C(1) - t1 * r_i_i_C(2) + t14 * t10 - t16 * t8, t8, t10 * t5, t3 * r_i_i_C(1) - t4 * r_i_i_C(2); t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t16 * t10 + t14 * t8, -t10, t8 * t5, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); 0, 0, -t6, -t11 * r_i_i_C(1) - t12 * r_i_i_C(2);];
	Ja_transl = t13;
else
	Ja_transl=NaN(3,4);
end