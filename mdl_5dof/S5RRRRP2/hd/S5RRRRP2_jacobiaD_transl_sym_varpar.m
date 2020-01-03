% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRP2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRP2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:12:14
	% EndTime: 2020-01-03 12:12:14
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:12:14
	% EndTime: 2020-01-03 12:12:14
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
	% StartTime: 2020-01-03 12:12:14
	% EndTime: 2020-01-03 12:12:14
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
	% StartTime: 2020-01-03 12:12:14
	% EndTime: 2020-01-03 12:12:14
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
	% StartTime: 2020-01-03 12:12:14
	% EndTime: 2020-01-03 12:12:14
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
	% StartTime: 2020-01-03 12:12:14
	% EndTime: 2020-01-03 12:12:14
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (264->42), mult. (202->56), div. (0->0), fcn. (131->8), ass. (0->35)
	t149 = pkin(4) + r_i_i_C(1);
	t120 = qJD(3) + qJD(4);
	t121 = qJD(1) + qJD(2);
	t122 = qJ(3) + qJ(4);
	t115 = sin(t122);
	t117 = cos(t122);
	t144 = r_i_i_C(2) * t117;
	t131 = -r_i_i_C(1) * t115 - t144;
	t124 = sin(qJ(3));
	t142 = pkin(3) * qJD(3);
	t133 = t124 * t142;
	t141 = t115 * t120;
	t148 = -pkin(4) * t141 - (-qJ(5) - pkin(8) - pkin(7)) * t121 + t120 * t131 - t133;
	t125 = cos(qJ(3));
	t147 = t125 * pkin(3) + t149 * t117 + pkin(2);
	t143 = pkin(1) * qJD(1);
	t140 = t117 * t120;
	t123 = qJ(1) + qJ(2);
	t118 = cos(t123);
	t139 = t118 * t120;
	t116 = sin(t123);
	t138 = t121 * t116;
	t137 = t121 * t118;
	t136 = r_i_i_C(1) * t140;
	t135 = r_i_i_C(2) * t141;
	t134 = r_i_i_C(2) * t138;
	t132 = t149 * t121;
	t130 = (-t124 * pkin(3) - pkin(4) * t115 + t131) * t121;
	t129 = (-t115 * t149 - t144) * t120;
	t127 = r_i_i_C(3) * t137 + t116 * qJD(5) + t115 * t134 + t148 * t118 - t147 * t138;
	t126 = r_i_i_C(3) * t138 + (-r_i_i_C(2) * t115 * t121 - qJD(5)) * t118 + t147 * t137 + t148 * t116;
	t107 = t118 * t136;
	t106 = t116 * t135;
	t105 = -pkin(4) * t140 - t125 * t142;
	t1 = [0, 0, t129 - t133, t129, 0; -sin(qJ(1)) * t143 + t127, t127, t106 + (t105 - t136) * t116 + t118 * t130, t106 - t118 * t115 * t132 + (-t116 * t120 * t149 - r_i_i_C(2) * t137) * t117, t138; cos(qJ(1)) * t143 + t126, t126, t107 + (-t105 - t135) * t118 + t116 * t130, t107 + (pkin(4) * t139 - t134) * t117 + (-r_i_i_C(2) * t139 - t116 * t132) * t115, -t137;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end