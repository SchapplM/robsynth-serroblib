% Analytische Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% Ja_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RRRPP5_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_jacobia_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPP5_jacobia_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPP5_jacobia_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_jacobia_transl_sym_varpar: pkin has to be [7x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:59:08
	% EndTime: 2019-12-31 20:59:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:59:08
	% EndTime: 2019-12-31 20:59:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (2->2), mult. (4->4), div. (0->0), fcn. (4->2), ass. (0->3)
	t2 = cos(qJ(1));
	t1 = sin(qJ(1));
	t3 = [-r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0, 0, 0, 0; r_i_i_C(1) * t2 - r_i_i_C(2) * t1, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_transl = t3;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:59:08
	% EndTime: 2019-12-31 20:59:08
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->6), mult. (24->10), div. (0->0), fcn. (24->4), ass. (0->9)
	t8 = pkin(6) + r_i_i_C(3);
	t1 = sin(qJ(2));
	t3 = cos(qJ(2));
	t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t1;
	t6 = -r_i_i_C(1) * t1 - r_i_i_C(2) * t3;
	t5 = pkin(1) + t7;
	t4 = cos(qJ(1));
	t2 = sin(qJ(1));
	t9 = [-t2 * t5 + t4 * t8, t6 * t4, 0, 0, 0; t2 * t8 + t4 * t5, t6 * t2, 0, 0, 0; 0, t7, 0, 0, 0;];
	Ja_transl = t9;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:59:08
	% EndTime: 2019-12-31 20:59:08
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (37->9), mult. (41->14), div. (0->0), fcn. (41->6), ass. (0->12)
	t6 = qJ(2) + qJ(3);
	t3 = sin(t6);
	t4 = cos(t6);
	t14 = t4 * r_i_i_C(1) - t3 * r_i_i_C(2);
	t17 = t14 + cos(qJ(2)) * pkin(2);
	t15 = r_i_i_C(3) + pkin(7) + pkin(6);
	t13 = -r_i_i_C(1) * t3 - r_i_i_C(2) * t4;
	t12 = pkin(1) + t17;
	t11 = -sin(qJ(2)) * pkin(2) + t13;
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t12 * t8 + t15 * t9, t11 * t9, t13 * t9, 0, 0; t12 * t9 + t15 * t8, t11 * t8, t13 * t8, 0, 0; 0, t17, t14, 0, 0;];
	Ja_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:59:08
	% EndTime: 2019-12-31 20:59:08
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (73->15), mult. (71->18), div. (0->0), fcn. (74->6), ass. (0->17)
	t32 = r_i_i_C(3) + qJ(4);
	t13 = qJ(2) + qJ(3);
	t10 = sin(t13);
	t11 = cos(t13);
	t27 = pkin(3) + r_i_i_C(1);
	t19 = t32 * t10 + t27 * t11;
	t31 = cos(qJ(2)) * pkin(2) + t19;
	t30 = t32 * t11;
	t28 = pkin(1) + t31;
	t15 = sin(qJ(1));
	t26 = t30 * t15;
	t16 = cos(qJ(1));
	t25 = t30 * t16;
	t23 = r_i_i_C(2) + pkin(7) + pkin(6);
	t20 = t27 * t10;
	t18 = -sin(qJ(2)) * pkin(2) - t20;
	t1 = [-t28 * t15 + t23 * t16, t18 * t16 + t25, -t16 * t20 + t25, t16 * t10, 0; t23 * t15 + t28 * t16, t18 * t15 + t26, -t15 * t20 + t26, t15 * t10, 0; 0, t31, t19, -t11, 0;];
	Ja_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 20:59:08
	% EndTime: 2019-12-31 20:59:08
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (92->16), mult. (87->18), div. (0->0), fcn. (92->6), ass. (0->17)
	t34 = r_i_i_C(2) + qJ(4);
	t14 = qJ(2) + qJ(3);
	t11 = sin(t14);
	t12 = cos(t14);
	t23 = pkin(3) + pkin(4) + r_i_i_C(1);
	t20 = t34 * t11 + t23 * t12;
	t33 = cos(qJ(2)) * pkin(2) + t20;
	t31 = t34 * t12;
	t29 = pkin(1) + t33;
	t16 = sin(qJ(1));
	t28 = t31 * t16;
	t17 = cos(qJ(1));
	t27 = t31 * t17;
	t22 = -r_i_i_C(3) - qJ(5) + pkin(7) + pkin(6);
	t21 = t23 * t11;
	t19 = -sin(qJ(2)) * pkin(2) - t21;
	t1 = [-t29 * t16 + t22 * t17, t19 * t17 + t27, -t17 * t21 + t27, t17 * t11, -t16; t22 * t16 + t29 * t17, t19 * t16 + t28, -t16 * t21 + t28, t16 * t11, t17; 0, t33, t20, -t12, 0;];
	Ja_transl = t1;
else
	Ja_transl=NaN(3,5);
end