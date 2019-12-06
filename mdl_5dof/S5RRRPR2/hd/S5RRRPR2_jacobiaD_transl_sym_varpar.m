% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:41:27
	% EndTime: 2019-12-05 18:41:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:41:27
	% EndTime: 2019-12-05 18:41:27
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
	% StartTime: 2019-12-05 18:41:27
	% EndTime: 2019-12-05 18:41:27
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
	% StartTime: 2019-12-05 18:41:27
	% EndTime: 2019-12-05 18:41:27
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (68->10), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->13)
	t24 = qJD(1) + qJD(2);
	t33 = pkin(2) * t24;
	t21 = qJD(3) + t24;
	t25 = qJ(1) + qJ(2);
	t23 = qJ(3) + t25;
	t32 = sin(t23) * t21;
	t31 = cos(t23) * t21;
	t30 = r_i_i_C(1) * t32 + r_i_i_C(2) * t31;
	t29 = pkin(1) * qJD(1);
	t28 = sin(t25) * t33 + t30;
	t27 = -r_i_i_C(1) * t31 + r_i_i_C(2) * t32;
	t26 = -cos(t25) * t33 + t27;
	t1 = [0, 0, 0, 0, 0; sin(qJ(1)) * t29 + t28, t28, t30, 0, 0; -cos(qJ(1)) * t29 + t26, t26, t27, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:41:27
	% EndTime: 2019-12-05 18:41:27
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (110->13), mult. (48->15), div. (0->0), fcn. (24->8), ass. (0->14)
	t30 = qJD(1) + qJD(2);
	t38 = pkin(2) * t30;
	t31 = qJ(1) + qJ(2);
	t29 = qJ(3) + t31;
	t25 = pkin(9) + t29;
	t27 = qJD(3) + t30;
	t37 = sin(t25) * t27;
	t36 = pkin(1) * qJD(1);
	t24 = cos(t25);
	t35 = r_i_i_C(1) * t37 + (sin(t29) * pkin(3) + r_i_i_C(2) * t24) * t27;
	t34 = sin(t31) * t38 + t35;
	t33 = r_i_i_C(2) * t37 + (-pkin(3) * cos(t29) - r_i_i_C(1) * t24) * t27;
	t32 = -cos(t31) * t38 + t33;
	t1 = [0, 0, 0, 0, 0; sin(qJ(1)) * t36 + t34, t34, t35, 0, 0; -cos(qJ(1)) * t36 + t32, t32, t33, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:41:27
	% EndTime: 2019-12-05 18:41:27
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (305->27), mult. (160->39), div. (0->0), fcn. (96->10), ass. (0->25)
	t105 = cos(qJ(5));
	t116 = qJD(5) * t105;
	t104 = sin(qJ(5));
	t102 = qJD(1) + qJD(2);
	t99 = qJD(3) + t102;
	t119 = t104 * t99;
	t103 = qJ(1) + qJ(2);
	t101 = qJ(3) + t103;
	t97 = pkin(9) + t101;
	t95 = sin(t97);
	t96 = cos(t97);
	t124 = t95 * t116 + t96 * t119;
	t117 = qJD(5) * t104;
	t118 = t105 * t99;
	t123 = t96 * t117 + t95 * t118;
	t122 = -pkin(8) - r_i_i_C(3);
	t121 = pkin(2) * t102;
	t120 = pkin(1) * qJD(1);
	t113 = t95 * t117;
	t110 = t96 * t116;
	t109 = r_i_i_C(2) * t110 + t123 * r_i_i_C(1) + (pkin(3) * sin(t101) + t122 * t96 + (-r_i_i_C(2) * t104 + pkin(4)) * t95) * t99;
	t108 = sin(t103) * t121 + t109;
	t107 = r_i_i_C(1) * t113 + (-pkin(3) * cos(t101) + (-r_i_i_C(1) * t105 - pkin(4)) * t96 + t122 * t95) * t99 + t124 * r_i_i_C(2);
	t106 = -cos(t103) * t121 + t107;
	t1 = [0, 0, 0, 0, (-r_i_i_C(1) * t104 - r_i_i_C(2) * t105) * qJD(5); sin(qJ(1)) * t120 + t108, t108, t109, 0, (t96 * t118 - t113) * r_i_i_C(2) + t124 * r_i_i_C(1); -cos(qJ(1)) * t120 + t106, t106, t107, 0, t123 * r_i_i_C(2) + (t95 * t119 - t110) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end