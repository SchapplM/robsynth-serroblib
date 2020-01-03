% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:59:31
	% EndTime: 2020-01-03 11:59:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:59:31
	% EndTime: 2020-01-03 11:59:31
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
	% StartTime: 2020-01-03 11:59:31
	% EndTime: 2020-01-03 11:59:31
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
	% StartTime: 2020-01-03 11:59:31
	% EndTime: 2020-01-03 11:59:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (42->9), mult. (28->11), div. (0->0), fcn. (14->6), ass. (0->9)
	t25 = pkin(1) * qJD(1);
	t22 = qJ(1) + qJ(2);
	t19 = pkin(8) + t22;
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
	% StartTime: 2020-01-03 11:59:31
	% EndTime: 2020-01-03 11:59:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (141->22), mult. (112->36), div. (0->0), fcn. (68->8), ass. (0->21)
	t113 = pkin(7) + r_i_i_C(3);
	t96 = sin(qJ(4));
	t106 = qJD(4) * t96;
	t94 = qJD(1) + qJD(2);
	t97 = cos(qJ(4));
	t108 = t94 * t97;
	t95 = qJ(1) + qJ(2);
	t92 = pkin(8) + t95;
	t90 = sin(t92);
	t91 = cos(t92);
	t112 = t90 * t106 - t91 * t108;
	t110 = t91 * t94;
	t109 = t94 * t96;
	t107 = pkin(1) * qJD(1);
	t105 = qJD(4) * t97;
	t104 = t90 * t109;
	t101 = (-r_i_i_C(1) * t96 - r_i_i_C(2) * t97) * qJD(4);
	t100 = -t90 * t105 - t91 * t109;
	t99 = pkin(3) * t110 - t112 * r_i_i_C(1) + t100 * r_i_i_C(2) + (cos(t95) * pkin(2) + t113 * t90) * t94;
	t98 = r_i_i_C(2) * t104 + t91 * t101 + (-pkin(2) * sin(t95) + (-r_i_i_C(1) * t97 - pkin(3)) * t90) * t94 + t113 * t110;
	t1 = [0, 0, 0, t101, 0; -sin(qJ(1)) * t107 + t98, t98, 0, t100 * r_i_i_C(1) + t112 * r_i_i_C(2), 0; cos(qJ(1)) * t107 + t99, t99, 0, (-t91 * t106 - t90 * t108) * r_i_i_C(2) + (t91 * t105 - t104) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:59:31
	% EndTime: 2020-01-03 11:59:31
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (195->24), mult. (148->31), div. (0->0), fcn. (95->8), ass. (0->23)
	t105 = cos(qJ(4));
	t117 = pkin(4) + r_i_i_C(1);
	t122 = t117 * t105;
	t121 = pkin(3) + t122;
	t104 = sin(qJ(4));
	t115 = r_i_i_C(2) * t104;
	t120 = qJD(4) * (-t115 + t122);
	t101 = qJD(1) + qJD(2);
	t119 = -t101 * t115 - qJD(5);
	t114 = pkin(1) * qJD(1);
	t102 = qJ(1) + qJ(2);
	t99 = pkin(8) + t102;
	t96 = sin(t99);
	t113 = t101 * t96;
	t97 = cos(t99);
	t112 = t101 * t97;
	t110 = -r_i_i_C(2) * t105 - t117 * t104;
	t109 = t101 * t110;
	t108 = t110 * qJD(4);
	t103 = -qJ(5) - pkin(7);
	t107 = r_i_i_C(3) * t113 + pkin(2) * t101 * cos(t102) + t119 * t97 + (-t101 * t103 + t108) * t96 + t121 * t112;
	t106 = r_i_i_C(3) * t112 + t97 * t108 + (-pkin(2) * sin(t102) - t103 * t97) * t101 + (-t121 * t101 - t119) * t96;
	t1 = [0, 0, 0, t108, 0; -sin(qJ(1)) * t114 + t106, t106, 0, t97 * t109 - t96 * t120, t113; cos(qJ(1)) * t114 + t107, t107, 0, t96 * t109 + t97 * t120, -t112;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end