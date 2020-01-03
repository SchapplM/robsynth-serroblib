% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR5
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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
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
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
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
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (71->13), mult. (58->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t76 = pkin(2) + r_i_i_C(1) * cos(pkin(9));
	t63 = qJD(1) + qJD(2);
	t75 = qJD(3) + r_i_i_C(2) * t63 * sin(pkin(9));
	t64 = qJ(1) + qJ(2);
	t61 = sin(t64);
	t73 = t63 * t61;
	t62 = cos(t64);
	t72 = t63 * t62;
	t71 = pkin(1) * qJD(1);
	t70 = qJ(3) * t63;
	t68 = r_i_i_C(3) * t72 + t75 * t61 + t62 * t70 - t76 * t73;
	t67 = r_i_i_C(3) * t73 + t61 * t70 - t75 * t62 + t76 * t72;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t71 + t68, t68, t73, 0, 0; cos(qJ(1)) * t71 + t67, t67, -t72, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (140->26), mult. (114->37), div. (0->0), fcn. (74->7), ass. (0->20)
	t98 = qJ(1) + qJ(2);
	t94 = sin(t98);
	t97 = qJD(1) + qJD(2);
	t110 = t97 * t94;
	t95 = cos(t98);
	t109 = t97 * t95;
	t108 = pkin(1) * qJD(1);
	t107 = qJD(4) * t94;
	t106 = qJD(4) * t95;
	t96 = pkin(9) + qJ(4);
	t92 = sin(t96);
	t105 = t92 * t110;
	t93 = cos(t96);
	t104 = t93 * t109;
	t103 = (-r_i_i_C(1) * t92 - r_i_i_C(2) * t93) * qJD(4);
	t102 = -t97 * (-pkin(7) - qJ(3)) + t103;
	t91 = cos(pkin(9)) * pkin(3) + pkin(2);
	t101 = r_i_i_C(2) * t105 + r_i_i_C(3) * t109 + t94 * qJD(3) + (-r_i_i_C(1) * t93 - t91) * t110 + t102 * t95;
	t100 = t91 * t109 + r_i_i_C(1) * t104 + r_i_i_C(3) * t110 + (-r_i_i_C(2) * t92 * t97 - qJD(3)) * t95 + t102 * t94;
	t1 = [0, 0, 0, t103, 0; -sin(qJ(1)) * t108 + t101, t101, t110, (t92 * t107 - t104) * r_i_i_C(2) + (-t93 * t107 - t92 * t109) * r_i_i_C(1), 0; cos(qJ(1)) * t108 + t100, t100, -t109, (-t92 * t106 - t93 * t110) * r_i_i_C(2) + (t93 * t106 - t105) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (265->38), mult. (172->49), div. (0->0), fcn. (113->9), ass. (0->33)
	t139 = pkin(1) * qJD(1);
	t138 = pkin(4) * qJD(4);
	t118 = pkin(9) + qJ(4);
	t114 = qJ(5) + t118;
	t110 = sin(t114);
	t119 = qJD(4) + qJD(5);
	t137 = t110 * t119;
	t111 = cos(t114);
	t136 = t111 * t119;
	t121 = qJ(1) + qJ(2);
	t115 = sin(t121);
	t120 = qJD(1) + qJD(2);
	t135 = t120 * t115;
	t116 = cos(t121);
	t134 = t120 * t116;
	t133 = r_i_i_C(1) * t136;
	t132 = r_i_i_C(2) * t137;
	t113 = cos(t118);
	t131 = t113 * t138;
	t130 = t110 * t135;
	t129 = t111 * t134;
	t128 = -r_i_i_C(1) * t110 - r_i_i_C(2) * t111;
	t127 = t128 * t119;
	t112 = sin(t118);
	t126 = (-pkin(4) * t112 + t128) * t120;
	t125 = -t112 * t138 + t127;
	t124 = -(-pkin(8) - pkin(7) - qJ(3)) * t120 + t125;
	t106 = pkin(4) * t113 + cos(pkin(9)) * pkin(3) + pkin(2);
	t123 = r_i_i_C(2) * t130 + r_i_i_C(3) * t134 + t115 * qJD(3) + (-r_i_i_C(1) * t111 - t106) * t135 + t124 * t116;
	t122 = t106 * t134 + r_i_i_C(1) * t129 + r_i_i_C(3) * t135 + (-r_i_i_C(2) * t110 * t120 - qJD(3)) * t116 + t124 * t115;
	t103 = t116 * t133;
	t102 = t115 * t132;
	t1 = [0, 0, 0, t125, t127; -sin(qJ(1)) * t139 + t123, t123, t135, t102 + (-t131 - t133) * t115 + t116 * t126, -r_i_i_C(2) * t129 + t102 + (-t110 * t134 - t115 * t136) * r_i_i_C(1); cos(qJ(1)) * t139 + t122, t122, -t134, t103 + (t131 - t132) * t116 + t115 * t126, -r_i_i_C(1) * t130 + t103 + (-t111 * t135 - t116 * t137) * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end