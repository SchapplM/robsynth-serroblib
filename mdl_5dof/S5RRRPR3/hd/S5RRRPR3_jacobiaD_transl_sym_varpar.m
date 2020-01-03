% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRPR3
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
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
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
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
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
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
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:14
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (159->22), mult. (140->28), div. (0->0), fcn. (91->8), ass. (0->22)
	t100 = qJ(3) + pkin(9);
	t96 = cos(t100);
	t122 = cos(qJ(3)) * pkin(3) + r_i_i_C(1) * t96;
	t121 = pkin(2) + t122;
	t95 = sin(t100);
	t117 = r_i_i_C(2) * t95;
	t99 = qJD(1) + qJD(2);
	t120 = t99 * t117 + qJD(4);
	t119 = qJD(3) * (-t117 + t122);
	t101 = qJ(1) + qJ(2);
	t97 = sin(t101);
	t115 = t99 * t97;
	t98 = cos(t101);
	t114 = t99 * t98;
	t113 = pkin(1) * qJD(1);
	t111 = -sin(qJ(3)) * pkin(3) - r_i_i_C(1) * t95 - r_i_i_C(2) * t96;
	t109 = t111 * t99;
	t108 = t111 * qJD(3);
	t107 = -(-qJ(4) - pkin(7)) * t99 + t108;
	t106 = r_i_i_C(3) * t114 + t107 * t98 - t121 * t115 + t120 * t97;
	t105 = r_i_i_C(3) * t115 + t107 * t97 + t121 * t114 - t120 * t98;
	t1 = [0, 0, t108, 0, 0; -sin(qJ(1)) * t113 + t106, t106, t98 * t109 - t119 * t97, t115, 0; cos(qJ(1)) * t113 + t105, t105, t97 * t109 + t119 * t98, -t114, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (274->39), mult. (188->49), div. (0->0), fcn. (122->10), ass. (0->32)
	t126 = qJ(3) + pkin(9);
	t112 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t126);
	t124 = qJD(3) + qJD(5);
	t120 = qJ(5) + t126;
	t116 = sin(t120);
	t117 = cos(t120);
	t135 = -r_i_i_C(1) * t116 - r_i_i_C(2) * t117;
	t134 = t135 * t124;
	t147 = t112 * qJD(3) + t134;
	t125 = qJD(1) + qJD(2);
	t146 = -(-pkin(8) - qJ(4) - pkin(7)) * t125 + t147;
	t145 = pkin(1) * qJD(1);
	t144 = t116 * t124;
	t143 = t117 * t124;
	t127 = qJ(1) + qJ(2);
	t121 = sin(t127);
	t142 = t125 * t121;
	t122 = cos(t127);
	t141 = t125 * t122;
	t140 = r_i_i_C(1) * t143;
	t139 = r_i_i_C(2) * t144;
	t138 = t116 * t142;
	t137 = t117 * t141;
	t136 = -cos(qJ(3)) * pkin(3) - pkin(4) * cos(t126);
	t133 = (t112 + t135) * t125;
	t111 = pkin(2) - t136;
	t131 = r_i_i_C(2) * t138 + r_i_i_C(3) * t141 + t121 * qJD(4) + (-r_i_i_C(1) * t117 - t111) * t142 + t146 * t122;
	t130 = t111 * t141 + r_i_i_C(1) * t137 + r_i_i_C(3) * t142 + (-r_i_i_C(2) * t116 * t125 - qJD(4)) * t122 + t146 * t121;
	t110 = t136 * qJD(3);
	t106 = t122 * t140;
	t105 = t121 * t139;
	t1 = [0, 0, t147, 0, t134; -sin(qJ(1)) * t145 + t131, t131, t105 + (t110 - t140) * t121 + t122 * t133, t142, -r_i_i_C(2) * t137 + t105 + (-t116 * t141 - t121 * t143) * r_i_i_C(1); cos(qJ(1)) * t145 + t130, t130, t106 + (-t110 - t139) * t122 + t121 * t133, -t141, -r_i_i_C(1) * t138 + t106 + (-t117 * t142 - t122 * t144) * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end