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
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
	% StartTime: 2019-12-05 18:31:00
	% EndTime: 2019-12-05 18:31:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:31:00
	% EndTime: 2019-12-05 18:31:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (r_i_i_C(1) * t5 + r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t6 + r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:31:00
	% EndTime: 2019-12-05 18:31:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t15 = qJD(1) + qJD(2);
	t16 = qJ(1) + qJ(2);
	t21 = sin(t16) * t15;
	t20 = cos(t16) * t15;
	t19 = r_i_i_C(1) * t21 + r_i_i_C(2) * t20;
	t18 = pkin(1) * qJD(1);
	t17 = -r_i_i_C(1) * t20 + r_i_i_C(2) * t21;
	t1 = [0, 0, 0, 0, 0; sin(qJ(1)) * t18 + t19, t19, 0, 0, 0; -cos(qJ(1)) * t18 + t17, t17, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:31:00
	% EndTime: 2019-12-05 18:31:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (42->9), mult. (28->12), div. (0->0), fcn. (14->6), ass. (0->9)
	t22 = qJ(1) + qJ(2);
	t19 = pkin(9) + t22;
	t21 = qJD(1) + qJD(2);
	t26 = sin(t19) * t21;
	t25 = pkin(1) * qJD(1);
	t18 = cos(t19);
	t24 = r_i_i_C(1) * t26 + (sin(t22) * pkin(2) + r_i_i_C(2) * t18) * t21;
	t23 = r_i_i_C(2) * t26 + (-pkin(2) * cos(t22) - r_i_i_C(1) * t18) * t21;
	t1 = [0, 0, 0, 0, 0; sin(qJ(1)) * t25 + t24, t24, 0, 0, 0; -cos(qJ(1)) * t25 + t23, t23, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:31:00
	% EndTime: 2019-12-05 18:31:00
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (96->13), mult. (44->15), div. (0->0), fcn. (22->8), ass. (0->13)
	t31 = qJ(1) + qJ(2);
	t28 = pkin(9) + t31;
	t26 = qJ(4) + t28;
	t30 = qJD(1) + qJD(2);
	t27 = qJD(4) + t30;
	t38 = sin(t26) * t27;
	t37 = cos(t26) * t27;
	t36 = r_i_i_C(1) * t38 + r_i_i_C(2) * t37;
	t35 = pkin(1) * qJD(1);
	t34 = t36 + (sin(t28) * pkin(3) + sin(t31) * pkin(2)) * t30;
	t33 = -r_i_i_C(1) * t37 + r_i_i_C(2) * t38;
	t32 = (-pkin(2) * cos(t31) - pkin(3) * cos(t28)) * t30 + t33;
	t1 = [0, 0, 0, 0, 0; sin(qJ(1)) * t35 + t34, t34, 0, t36, 0; -cos(qJ(1)) * t35 + t32, t32, 0, t33, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:31:00
	% EndTime: 2019-12-05 18:31:00
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (291->27), mult. (156->40), div. (0->0), fcn. (94->10), ass. (0->24)
	t105 = cos(qJ(5));
	t116 = qJD(5) * t105;
	t104 = sin(qJ(5));
	t102 = qJD(1) + qJD(2);
	t99 = qJD(4) + t102;
	t119 = t104 * t99;
	t103 = qJ(1) + qJ(2);
	t100 = pkin(9) + t103;
	t98 = qJ(4) + t100;
	t95 = sin(t98);
	t96 = cos(t98);
	t123 = t95 * t116 + t96 * t119;
	t117 = qJD(5) * t104;
	t118 = t105 * t99;
	t122 = t96 * t117 + t95 * t118;
	t121 = -pkin(8) - r_i_i_C(3);
	t120 = pkin(1) * qJD(1);
	t113 = t95 * t117;
	t110 = t96 * t116;
	t109 = r_i_i_C(2) * t110 + (t121 * t96 + (-r_i_i_C(2) * t104 + pkin(4)) * t95) * t99 + t122 * r_i_i_C(1);
	t108 = r_i_i_C(1) * t113 + ((-r_i_i_C(1) * t105 - pkin(4)) * t96 + t121 * t95) * t99 + t123 * r_i_i_C(2);
	t107 = t109 + (sin(t100) * pkin(3) + sin(t103) * pkin(2)) * t102;
	t106 = (-pkin(2) * cos(t103) - pkin(3) * cos(t100)) * t102 + t108;
	t1 = [0, 0, 0, 0, (-r_i_i_C(1) * t104 - r_i_i_C(2) * t105) * qJD(5); sin(qJ(1)) * t120 + t107, t107, 0, t109, (t96 * t118 - t113) * r_i_i_C(2) + t123 * r_i_i_C(1); -cos(qJ(1)) * t120 + t106, t106, 0, t108, t122 * r_i_i_C(2) + (t95 * t119 - t110) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end