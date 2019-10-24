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
% Datum: 2019-10-24 10:48
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
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
	% StartTime: 2019-10-24 10:48:01
	% EndTime: 2019-10-24 10:48:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:48:01
	% EndTime: 2019-10-24 10:48:01
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
	% StartTime: 2019-10-24 10:48:01
	% EndTime: 2019-10-24 10:48:01
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
	% StartTime: 2019-10-24 10:48:01
	% EndTime: 2019-10-24 10:48:01
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (42->9), mult. (28->12), div. (0->0), fcn. (14->6), ass. (0->9)
	t22 = qJ(1) + qJ(2);
	t19 = pkin(8) + t22;
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
	% StartTime: 2019-10-24 10:48:01
	% EndTime: 2019-10-24 10:48:02
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (141->23), mult. (112->36), div. (0->0), fcn. (68->8), ass. (0->20)
	t95 = sin(qJ(4));
	t106 = qJD(4) * t95;
	t93 = qJD(1) + qJD(2);
	t96 = cos(qJ(4));
	t108 = t93 * t96;
	t94 = qJ(1) + qJ(2);
	t91 = pkin(8) + t94;
	t89 = sin(t91);
	t90 = cos(t91);
	t112 = t90 * t106 + t89 * t108;
	t105 = qJD(4) * t96;
	t109 = t93 * t95;
	t111 = t89 * t105 + t90 * t109;
	t110 = -pkin(7) - r_i_i_C(3);
	t107 = pkin(1) * qJD(1);
	t102 = t89 * t106;
	t99 = t90 * t105;
	t98 = r_i_i_C(2) * t99 + t112 * r_i_i_C(1) + (pkin(2) * sin(t94) + t110 * t90 + (-r_i_i_C(2) * t95 + pkin(3)) * t89) * t93;
	t97 = r_i_i_C(1) * t102 + (-pkin(2) * cos(t94) + (-r_i_i_C(1) * t96 - pkin(3)) * t90 + t110 * t89) * t93 + t111 * r_i_i_C(2);
	t1 = [0, 0, 0, (-r_i_i_C(1) * t95 - r_i_i_C(2) * t96) * qJD(4), 0; sin(qJ(1)) * t107 + t98, t98, 0, (t90 * t108 - t102) * r_i_i_C(2) + t111 * r_i_i_C(1), 0; -cos(qJ(1)) * t107 + t97, t97, 0, t112 * r_i_i_C(2) + (t89 * t109 - t99) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:48:01
	% EndTime: 2019-10-24 10:48:01
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (195->24), mult. (148->29), div. (0->0), fcn. (95->8), ass. (0->24)
	t109 = sin(qJ(4));
	t110 = cos(qJ(4));
	t126 = pkin(4) + r_i_i_C(1);
	t128 = r_i_i_C(2) * t110 + t126 * t109;
	t138 = t128 * qJD(4);
	t137 = t126 * t110;
	t106 = qJD(1) + qJD(2);
	t135 = t106 * (-qJ(5) - pkin(7)) + t138;
	t134 = pkin(2) * t106;
	t133 = pkin(3) + t137;
	t124 = r_i_i_C(2) * t109;
	t131 = qJD(4) * (-t124 + t137);
	t130 = -t106 * t124 - qJD(5);
	t122 = pkin(1) * qJD(1);
	t107 = qJ(1) + qJ(2);
	t104 = pkin(8) + t107;
	t101 = sin(t104);
	t121 = t106 * t101;
	t102 = cos(t104);
	t120 = t106 * t102;
	t113 = t128 * t106;
	t112 = sin(t107) * t134 - r_i_i_C(3) * t120 + t130 * t101 + t133 * t121 + t135 * t102;
	t111 = -cos(t107) * t134 + (-t133 * t106 - t130) * t102 + (-r_i_i_C(3) * t106 + t135) * t101;
	t1 = [0, 0, 0, -t138, 0; sin(qJ(1)) * t122 + t112, t112, 0, t101 * t131 + t102 * t113, -t121; -cos(qJ(1)) * t122 + t111, t111, 0, t101 * t113 - t102 * t131, t120;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end