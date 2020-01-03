% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:16:04
	% EndTime: 2020-01-03 12:16:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:16:04
	% EndTime: 2020-01-03 12:16:04
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
	% StartTime: 2020-01-03 12:16:04
	% EndTime: 2020-01-03 12:16:04
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
	% StartTime: 2020-01-03 12:16:04
	% EndTime: 2020-01-03 12:16:04
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (93->18), mult. (104->32), div. (0->0), fcn. (64->6), ass. (0->21)
	t109 = pkin(7) + r_i_i_C(3);
	t91 = cos(qJ(3));
	t100 = qJD(3) * t91;
	t88 = qJD(1) + qJD(2);
	t90 = sin(qJ(3));
	t104 = t88 * t90;
	t89 = qJ(1) + qJ(2);
	t86 = sin(t89);
	t87 = cos(t89);
	t108 = t87 * t100 - t86 * t104;
	t101 = qJD(3) * t90;
	t103 = t88 * t91;
	t107 = t86 * t101 - t87 * t103;
	t106 = t86 * t88;
	t105 = t87 * t88;
	t102 = pkin(1) * qJD(1);
	t95 = -t87 * t101 - t86 * t103;
	t94 = -t86 * t100 - t87 * t104;
	t93 = pkin(2) * t105 - t107 * r_i_i_C(1) + t94 * r_i_i_C(2) + t109 * t106;
	t92 = -pkin(2) * t106 + t95 * r_i_i_C(1) - t108 * r_i_i_C(2) + t109 * t105;
	t1 = [0, 0, (-r_i_i_C(1) * t90 - r_i_i_C(2) * t91) * qJD(3), 0, 0; -sin(qJ(1)) * t102 + t92, t92, t94 * r_i_i_C(1) + t107 * r_i_i_C(2), 0, 0; cos(qJ(1)) * t102 + t93, t93, t108 * r_i_i_C(1) + t95 * r_i_i_C(2), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:16:04
	% EndTime: 2020-01-03 12:16:04
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (195->32), mult. (162->45), div. (0->0), fcn. (103->8), ass. (0->33)
	t137 = pkin(1) * qJD(1);
	t136 = pkin(3) * qJD(3);
	t114 = qJ(3) + qJ(4);
	t108 = sin(t114);
	t112 = qJD(3) + qJD(4);
	t135 = t108 * t112;
	t115 = qJ(1) + qJ(2);
	t109 = sin(t115);
	t113 = qJD(1) + qJD(2);
	t134 = t109 * t113;
	t110 = cos(t114);
	t133 = t110 * t112;
	t111 = cos(t115);
	t132 = t111 * t113;
	t131 = r_i_i_C(1) * t133;
	t130 = r_i_i_C(2) * t135;
	t117 = cos(qJ(3));
	t129 = t117 * t136;
	t128 = t108 * t134;
	t127 = t108 * t132;
	t126 = t110 * t132;
	t125 = -r_i_i_C(1) * t108 - r_i_i_C(2) * t110;
	t124 = t125 * t112;
	t116 = sin(qJ(3));
	t123 = (-pkin(3) * t116 + t125) * t113;
	t122 = -t116 * t136 + t124;
	t121 = -t113 * (-pkin(8) - pkin(7)) + t122;
	t107 = t117 * pkin(3) + pkin(2);
	t120 = r_i_i_C(1) * t126 - r_i_i_C(2) * t127 + r_i_i_C(3) * t134 + t107 * t132 + t121 * t109;
	t119 = r_i_i_C(2) * t128 + r_i_i_C(3) * t132 + (-r_i_i_C(1) * t110 - t107) * t134 + t121 * t111;
	t102 = t111 * t131;
	t101 = t109 * t130;
	t1 = [0, 0, t122, t124, 0; -sin(qJ(1)) * t137 + t119, t119, t101 + (-t129 - t131) * t109 + t111 * t123, -r_i_i_C(2) * t126 + t101 + (-t109 * t133 - t127) * r_i_i_C(1), 0; cos(qJ(1)) * t137 + t120, t120, t102 + (t129 - t130) * t111 + t109 * t123, -r_i_i_C(1) * t128 + t102 + (-t110 * t134 - t111 * t135) * r_i_i_C(2), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:16:04
	% EndTime: 2020-01-03 12:16:04
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (358->45), mult. (224->55), div. (0->0), fcn. (143->10), ass. (0->42)
	t130 = sin(qJ(3));
	t128 = qJ(3) + qJ(4);
	t120 = sin(t128);
	t125 = qJD(3) + qJD(4);
	t119 = qJD(5) + t125;
	t124 = qJ(5) + t128;
	t117 = sin(t124);
	t118 = cos(t124);
	t139 = -r_i_i_C(1) * t117 - r_i_i_C(2) * t118;
	t138 = t139 * t119;
	t154 = pkin(4) * t125;
	t134 = -t120 * t154 + t138;
	t152 = pkin(3) * qJD(3);
	t158 = -t130 * t152 + t134;
	t126 = qJD(1) + qJD(2);
	t157 = -t126 * (-pkin(9) - pkin(8) - pkin(7)) + t158;
	t156 = -pkin(4) * t120 + t139;
	t153 = pkin(1) * qJD(1);
	t151 = t117 * t119;
	t150 = t118 * t119;
	t129 = qJ(1) + qJ(2);
	t121 = sin(t129);
	t149 = t121 * t126;
	t123 = cos(t129);
	t148 = t123 * t126;
	t122 = cos(t128);
	t146 = t122 * t154;
	t145 = r_i_i_C(1) * t150;
	t144 = r_i_i_C(2) * t151;
	t142 = t117 * t149;
	t141 = t117 * t148;
	t140 = t118 * t148;
	t137 = (-t130 * pkin(3) + t156) * t126;
	t136 = t156 * t126;
	t131 = cos(qJ(3));
	t113 = t131 * pkin(3) + pkin(4) * t122 + pkin(2);
	t133 = r_i_i_C(1) * t140 - r_i_i_C(2) * t141 + r_i_i_C(3) * t149 + t113 * t148 + t157 * t121;
	t132 = r_i_i_C(2) * t142 + r_i_i_C(3) * t148 + (-r_i_i_C(1) * t118 - t113) * t149 + t157 * t123;
	t110 = -t131 * t152 - t146;
	t108 = t123 * t145;
	t107 = t121 * t144;
	t1 = [0, 0, t158, t134, t138; -sin(qJ(1)) * t153 + t132, t132, t107 + (t110 - t145) * t121 + t123 * t137, t107 + (-t145 - t146) * t121 + t123 * t136, -r_i_i_C(2) * t140 + t107 + (-t121 * t150 - t141) * r_i_i_C(1); cos(qJ(1)) * t153 + t133, t133, t108 + (-t110 - t144) * t123 + t121 * t137, t108 + (-t144 + t146) * t123 + t121 * t136, -r_i_i_C(1) * t142 + t108 + (-t118 * t149 - t123 * t151) * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end