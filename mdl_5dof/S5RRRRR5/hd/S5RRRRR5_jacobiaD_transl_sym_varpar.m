% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRRRR5
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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
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
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
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
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (68->10), mult. (36->12), div. (0->0), fcn. (18->6), ass. (0->13)
	t24 = qJD(1) + qJD(2);
	t31 = pkin(2) * t24;
	t30 = pkin(1) * qJD(1);
	t25 = qJ(1) + qJ(2);
	t23 = qJ(3) + t25;
	t19 = sin(t23);
	t20 = cos(t23);
	t21 = qJD(3) + t24;
	t29 = (r_i_i_C(1) * t20 - r_i_i_C(2) * t19) * t21;
	t28 = cos(t25) * t31 + t29;
	t27 = (-r_i_i_C(1) * t19 - r_i_i_C(2) * t20) * t21;
	t26 = -sin(t25) * t31 + t27;
	t1 = [0, 0, 0, 0, 0; -sin(qJ(1)) * t30 + t26, t26, t27, 0, 0; cos(qJ(1)) * t30 + t28, t28, t29, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (225->22), mult. (148->35), div. (0->0), fcn. (90->8), ass. (0->26)
	t121 = pkin(8) + r_i_i_C(3);
	t100 = cos(qJ(4));
	t111 = qJD(4) * t100;
	t97 = qJD(1) + qJD(2);
	t94 = qJD(3) + t97;
	t99 = sin(qJ(4));
	t115 = t94 * t99;
	t98 = qJ(1) + qJ(2);
	t96 = qJ(3) + t98;
	t92 = sin(t96);
	t93 = cos(t96);
	t120 = t93 * t111 - t92 * t115;
	t112 = qJD(4) * t99;
	t113 = t100 * t94;
	t119 = t92 * t112 - t93 * t113;
	t118 = pkin(2) * t97;
	t117 = t92 * t94;
	t116 = t93 * t94;
	t114 = pkin(1) * qJD(1);
	t106 = -t93 * t112 - t92 * t113;
	t105 = -t92 * t111 - t93 * t115;
	t104 = pkin(3) * t116 - t119 * r_i_i_C(1) + t105 * r_i_i_C(2) + t117 * t121;
	t103 = cos(t98) * t118 + t104;
	t102 = -pkin(3) * t117 + t106 * r_i_i_C(1) - t120 * r_i_i_C(2) + t116 * t121;
	t101 = -sin(t98) * t118 + t102;
	t1 = [0, 0, 0, (-r_i_i_C(1) * t99 - r_i_i_C(2) * t100) * qJD(4), 0; -sin(qJ(1)) * t114 + t101, t101, t102, t105 * r_i_i_C(1) + r_i_i_C(2) * t119, 0; cos(qJ(1)) * t114 + t103, t103, t104, r_i_i_C(1) * t120 + t106 * r_i_i_C(2), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:13:57
	% EndTime: 2020-01-03 12:13:57
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (371->36), mult. (214->50), div. (0->0), fcn. (135->10), ass. (0->40)
	t122 = qJD(1) + qJD(2);
	t151 = pkin(2) * t122;
	t150 = pkin(1) * qJD(1);
	t149 = pkin(4) * qJD(4);
	t124 = qJ(1) + qJ(2);
	t120 = qJ(3) + t124;
	t113 = sin(t120);
	t116 = qJD(3) + t122;
	t148 = t113 * t116;
	t114 = cos(t120);
	t147 = t114 * t116;
	t123 = qJ(4) + qJ(5);
	t117 = sin(t123);
	t146 = t116 * t117;
	t119 = cos(t123);
	t145 = t116 * t119;
	t121 = qJD(4) + qJD(5);
	t144 = t117 * t121;
	t143 = t119 * t121;
	t142 = r_i_i_C(1) * t143;
	t141 = r_i_i_C(2) * t144;
	t126 = cos(qJ(4));
	t140 = t126 * t149;
	t139 = t113 * t146;
	t138 = t114 * t146;
	t137 = t114 * t145;
	t136 = -r_i_i_C(1) * t117 - r_i_i_C(2) * t119;
	t135 = t136 * t121;
	t125 = sin(qJ(4));
	t134 = (-pkin(4) * t125 + t136) * t116;
	t133 = -t125 * t149 + t135;
	t132 = -t116 * (-pkin(9) - pkin(8)) + t133;
	t115 = t126 * pkin(4) + pkin(3);
	t131 = r_i_i_C(1) * t137 - r_i_i_C(2) * t138 + r_i_i_C(3) * t148 + t132 * t113 + t115 * t147;
	t130 = cos(t124) * t151 + t131;
	t129 = r_i_i_C(2) * t139 + r_i_i_C(3) * t147 + (-r_i_i_C(1) * t119 - t115) * t148 + t132 * t114;
	t128 = -sin(t124) * t151 + t129;
	t109 = t114 * t142;
	t108 = t113 * t141;
	t1 = [0, 0, 0, t133, t135; -sin(qJ(1)) * t150 + t128, t128, t129, t108 + (-t140 - t142) * t113 + t114 * t134, -r_i_i_C(2) * t137 + t108 + (-t113 * t143 - t138) * r_i_i_C(1); cos(qJ(1)) * t150 + t130, t130, t131, t109 + (t140 - t141) * t114 + t113 * t134, -r_i_i_C(1) * t139 + t109 + (-t113 * t145 - t114 * t144) * r_i_i_C(2);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end