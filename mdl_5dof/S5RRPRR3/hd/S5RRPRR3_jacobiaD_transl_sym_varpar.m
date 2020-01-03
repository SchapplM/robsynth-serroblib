% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR3
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:00:59
	% EndTime: 2020-01-03 12:00:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:00:59
	% EndTime: 2020-01-03 12:00:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:00:59
	% EndTime: 2020-01-03 12:00:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t19 = pkin(1) * qJD(1);
	t16 = qJ(1) + qJ(2);
	t13 = sin(t16);
	t14 = cos(t16);
	t15 = qJD(1) + qJD(2);
	t18 = (r_i_i_C(1) * t14 - r_i_i_C(2) * t13) * t15;
	t17 = (-r_i_i_C(1) * t13 - r_i_i_C(2) * t14) * t15;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t19 + t17, t17, 0, 0, 0; cos(qJ(1)) * t19 + t18, t18, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:00:59
	% EndTime: 2020-01-03 12:00:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (42->9), mult. (28->11), div. (0->0), fcn. (14->6), ass. (0->9)
	t25 = pkin(1) * qJD(1);
	t22 = qJ(1) + qJ(2);
	t19 = pkin(9) + t22;
	t17 = sin(t19);
	t18 = cos(t19);
	t21 = qJD(1) + qJD(2);
	t24 = (cos(t22) * pkin(2) + r_i_i_C(1) * t18 - r_i_i_C(2) * t17) * t21;
	t23 = (-pkin(2) * sin(t22) - r_i_i_C(1) * t17 - r_i_i_C(2) * t18) * t21;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t25 + t23, t23, 0, 0, 0; cos(qJ(1)) * t25 + t24, t24, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:00:59
	% EndTime: 2020-01-03 12:00:59
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (96->13), mult. (44->15), div. (0->0), fcn. (22->8), ass. (0->13)
	t36 = pkin(1) * qJD(1);
	t31 = qJ(1) + qJ(2);
	t30 = qJD(1) + qJD(2);
	t28 = pkin(9) + t31;
	t26 = qJ(4) + t28;
	t23 = sin(t26);
	t24 = cos(t26);
	t27 = qJD(4) + t30;
	t35 = (r_i_i_C(1) * t24 - r_i_i_C(2) * t23) * t27;
	t34 = t35 + (cos(t28) * pkin(3) + cos(t31) * pkin(2)) * t30;
	t33 = (-r_i_i_C(1) * t23 - r_i_i_C(2) * t24) * t27;
	t32 = (-pkin(2) * sin(t31) - pkin(3) * sin(t28)) * t30 + t33;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t36 + t32, t32, 0, t33, 0; cos(qJ(1)) * t36 + t34, t34, 0, t35, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:00:59
	% EndTime: 2020-01-03 12:00:59
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (291->25), mult. (156->38), div. (0->0), fcn. (94->10), ass. (0->26)
	t126 = pkin(8) + r_i_i_C(3);
	t106 = cos(qJ(5));
	t117 = qJD(5) * t106;
	t103 = qJD(1) + qJD(2);
	t100 = qJD(4) + t103;
	t105 = sin(qJ(5));
	t120 = t100 * t105;
	t104 = qJ(1) + qJ(2);
	t101 = pkin(9) + t104;
	t99 = qJ(4) + t101;
	t96 = sin(t99);
	t97 = cos(t99);
	t125 = t97 * t117 - t96 * t120;
	t118 = qJD(5) * t105;
	t119 = t100 * t106;
	t124 = t96 * t118 - t97 * t119;
	t123 = pkin(1) * qJD(1);
	t122 = t100 * t96;
	t121 = t100 * t97;
	t112 = -t96 * t117 - t97 * t120;
	t111 = -t97 * t118 - t96 * t119;
	t110 = pkin(4) * t121 - t124 * r_i_i_C(1) + t112 * r_i_i_C(2) + t126 * t122;
	t109 = t110 + (cos(t101) * pkin(3) + cos(t104) * pkin(2)) * t103;
	t108 = -pkin(4) * t122 + t111 * r_i_i_C(1) - t125 * r_i_i_C(2) + t126 * t121;
	t107 = (-pkin(2) * sin(t104) - pkin(3) * sin(t101)) * t103 + t108;
	t1 = [0, 0, 0, 0, (-r_i_i_C(1) * t105 - r_i_i_C(2) * t106) * qJD(5); -sin(qJ(1)) * t123 + t107, t107, 0, t108, t112 * r_i_i_C(1) + t124 * r_i_i_C(2); cos(qJ(1)) * t123 + t109, t109, 0, t110, t125 * r_i_i_C(1) + t111 * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end