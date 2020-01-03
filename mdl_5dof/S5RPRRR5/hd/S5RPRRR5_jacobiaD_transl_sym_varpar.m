% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t6 = cos(qJ(1));
	t5 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; (-r_i_i_C(1) * t5 - r_i_i_C(2) * t6) * qJD(1), 0, 0, 0, 0; (r_i_i_C(1) * t6 - r_i_i_C(2) * t5) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t9 = qJ(1) + pkin(9);
	t8 = cos(t9);
	t7 = sin(t9);
	t1 = [0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t7 - r_i_i_C(2) * t8) * qJD(1), 0, 0, 0, 0; (cos(qJ(1)) * pkin(1) + r_i_i_C(1) * t8 - r_i_i_C(2) * t7) * qJD(1), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (34->9), mult. (24->12), div. (0->0), fcn. (12->6), ass. (0->8)
	t18 = qJ(1) + pkin(9);
	t16 = qJ(3) + t18;
	t14 = sin(t16);
	t15 = cos(t16);
	t17 = qJD(1) + qJD(3);
	t20 = (r_i_i_C(1) * t15 - r_i_i_C(2) * t14) * t17;
	t19 = (-r_i_i_C(1) * t14 - r_i_i_C(2) * t15) * t17;
	t1 = [0, 0, 0, 0, 0; t19 + (-sin(t18) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1), 0, t19, 0, 0; (cos(t18) * pkin(2) + cos(qJ(1)) * pkin(1)) * qJD(1) + t20, 0, t20, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (133->21), mult. (108->35), div. (0->0), fcn. (66->8), ass. (0->21)
	t110 = pkin(7) + r_i_i_C(3);
	t93 = cos(qJ(4));
	t102 = qJD(4) * t93;
	t90 = qJD(1) + qJD(3);
	t92 = sin(qJ(4));
	t105 = t90 * t92;
	t91 = qJ(1) + pkin(9);
	t89 = qJ(3) + t91;
	t87 = sin(t89);
	t88 = cos(t89);
	t109 = t88 * t102 - t87 * t105;
	t103 = qJD(4) * t92;
	t104 = t90 * t93;
	t108 = t87 * t103 - t88 * t104;
	t107 = t87 * t90;
	t106 = t88 * t90;
	t97 = -t88 * t103 - t87 * t104;
	t96 = -t87 * t102 - t88 * t105;
	t95 = pkin(3) * t106 - t108 * r_i_i_C(1) + t96 * r_i_i_C(2) + t107 * t110;
	t94 = -pkin(3) * t107 + t97 * r_i_i_C(1) - t109 * r_i_i_C(2) + t106 * t110;
	t1 = [0, 0, 0, (-r_i_i_C(1) * t92 - r_i_i_C(2) * t93) * qJD(4), 0; (-sin(qJ(1)) * pkin(1) - sin(t91) * pkin(2)) * qJD(1) + t94, 0, t94, t96 * r_i_i_C(1) + r_i_i_C(2) * t108, 0; (cos(t91) * pkin(2) + cos(qJ(1)) * pkin(1)) * qJD(1) + t95, 0, t95, r_i_i_C(1) * t109 + t97 * r_i_i_C(2), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:54:46
	% EndTime: 2020-01-03 11:54:46
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (251->35), mult. (166->50), div. (0->0), fcn. (105->10), ass. (0->35)
	t140 = pkin(4) * qJD(4);
	t116 = qJ(1) + pkin(9);
	t111 = qJ(3) + t116;
	t108 = sin(t111);
	t115 = qJD(1) + qJD(3);
	t139 = t108 * t115;
	t109 = cos(t111);
	t138 = t109 * t115;
	t117 = qJ(4) + qJ(5);
	t112 = sin(t117);
	t114 = qJD(4) + qJD(5);
	t137 = t112 * t114;
	t136 = t112 * t115;
	t113 = cos(t117);
	t135 = t113 * t114;
	t134 = t113 * t115;
	t133 = r_i_i_C(1) * t135;
	t132 = r_i_i_C(2) * t137;
	t119 = cos(qJ(4));
	t131 = t119 * t140;
	t130 = t108 * t136;
	t129 = t109 * t136;
	t128 = t109 * t134;
	t127 = -r_i_i_C(1) * t112 - r_i_i_C(2) * t113;
	t126 = t127 * t114;
	t118 = sin(qJ(4));
	t125 = (-pkin(4) * t118 + t127) * t115;
	t124 = -t118 * t140 + t126;
	t123 = -t115 * (-pkin(8) - pkin(7)) + t124;
	t110 = t119 * pkin(4) + pkin(3);
	t122 = r_i_i_C(1) * t128 - r_i_i_C(2) * t129 + r_i_i_C(3) * t139 + t123 * t108 + t110 * t138;
	t121 = r_i_i_C(2) * t130 + r_i_i_C(3) * t138 + (-r_i_i_C(1) * t113 - t110) * t139 + t123 * t109;
	t103 = t109 * t133;
	t102 = t108 * t132;
	t1 = [0, 0, 0, t124, t126; (-sin(t116) * pkin(2) - sin(qJ(1)) * pkin(1)) * qJD(1) + t121, 0, t121, t102 + (-t131 - t133) * t108 + t109 * t125, -r_i_i_C(2) * t128 + t102 + (-t108 * t135 - t129) * r_i_i_C(1); (cos(t116) * pkin(2) + cos(qJ(1)) * pkin(1)) * qJD(1) + t122, 0, t122, t103 + (t131 - t132) * t109 + t108 * t125, -r_i_i_C(1) * t130 + t103 + (-t108 * t134 - t109 * t137) * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end