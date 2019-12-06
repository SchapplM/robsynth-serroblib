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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
	% StartTime: 2019-12-05 18:34:48
	% EndTime: 2019-12-05 18:34:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:48
	% EndTime: 2019-12-05 18:34:48
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
	% StartTime: 2019-12-05 18:34:49
	% EndTime: 2019-12-05 18:34:49
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
	% StartTime: 2019-12-05 18:34:49
	% EndTime: 2019-12-05 18:34:49
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (71->13), mult. (58->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t69 = pkin(2) + r_i_i_C(1) * cos(pkin(9));
	t67 = r_i_i_C(2) * sin(pkin(9));
	t58 = qJ(1) + qJ(2);
	t55 = sin(t58);
	t57 = qJD(1) + qJD(2);
	t66 = t57 * t55;
	t56 = cos(t58);
	t65 = t57 * t56;
	t64 = -r_i_i_C(3) - qJ(3);
	t63 = pkin(1) * qJD(1);
	t62 = t65 * t67 + t56 * qJD(3) + (t64 * t55 - t56 * t69) * t57;
	t61 = -t55 * qJD(3) + (-t55 * t67 + t56 * t64) * t57 + t69 * t66;
	t1 = [0, 0, 0, 0, 0; sin(qJ(1)) * t63 + t61, t61, -t66, 0, 0; -cos(qJ(1)) * t63 + t62, t62, t65, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:49
	% EndTime: 2019-12-05 18:34:49
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (140->27), mult. (114->41), div. (0->0), fcn. (74->7), ass. (0->21)
	t101 = qJ(1) + qJ(2);
	t97 = sin(t101);
	t113 = qJD(4) * t97;
	t100 = qJD(1) + qJD(2);
	t98 = cos(t101);
	t114 = t100 * t98;
	t99 = pkin(9) + qJ(4);
	t95 = sin(t99);
	t96 = cos(t99);
	t118 = t96 * t113 + t95 * t114;
	t112 = qJD(4) * t98;
	t115 = t100 * t97;
	t117 = t95 * t112 + t96 * t115;
	t116 = pkin(1) * qJD(1);
	t111 = t100 * (-pkin(7) - qJ(3));
	t108 = t95 * t113;
	t105 = t96 * t112;
	t94 = cos(pkin(9)) * pkin(3) + pkin(2);
	t104 = r_i_i_C(1) * t108 + t97 * t111 + t98 * qJD(3) + (-r_i_i_C(3) * t97 + (-r_i_i_C(1) * t96 - t94) * t98) * t100 + t118 * r_i_i_C(2);
	t103 = -t97 * qJD(3) + t94 * t115 + r_i_i_C(2) * t105 + t98 * t111 + (-r_i_i_C(2) * t95 * t97 - r_i_i_C(3) * t98) * t100 + t117 * r_i_i_C(1);
	t1 = [0, 0, 0, (-r_i_i_C(1) * t95 - r_i_i_C(2) * t96) * qJD(4), 0; sin(qJ(1)) * t116 + t103, t103, -t115, (t96 * t114 - t108) * r_i_i_C(2) + t118 * r_i_i_C(1), 0; -cos(qJ(1)) * t116 + t104, t104, t114, t117 * r_i_i_C(2) + (t95 * t115 - t105) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:49
	% EndTime: 2019-12-05 18:34:49
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (265->37), mult. (172->49), div. (0->0), fcn. (113->9), ass. (0->32)
	t127 = pkin(9) + qJ(4);
	t123 = qJ(5) + t127;
	t119 = sin(t123);
	t120 = cos(t123);
	t130 = qJ(1) + qJ(2);
	t125 = cos(t130);
	t129 = qJD(1) + qJD(2);
	t144 = t129 * t125;
	t124 = sin(t130);
	t128 = qJD(4) + qJD(5);
	t147 = t124 * t128;
	t152 = t119 * t144 + t120 * t147;
	t145 = t129 * t124;
	t146 = t125 * t128;
	t151 = t119 * t146 + t120 * t145;
	t150 = r_i_i_C(1) * t119;
	t149 = r_i_i_C(2) * t120;
	t148 = pkin(1) * qJD(1);
	t122 = cos(t127);
	t143 = qJD(4) * t122;
	t121 = sin(t127);
	t142 = pkin(4) * qJD(4) * t121;
	t141 = t119 * t147;
	t136 = t120 * t146;
	t135 = (-t149 - t150) * t128;
	t134 = -r_i_i_C(1) * t136 + t151 * r_i_i_C(2) + t145 * t150;
	t133 = t152 * r_i_i_C(1) - r_i_i_C(2) * t141 + t144 * t149;
	t113 = pkin(4) * t122 + cos(pkin(9)) * pkin(3) + pkin(2);
	t126 = -pkin(8) - pkin(7) - qJ(3);
	t132 = r_i_i_C(1) * t141 + t124 * t142 + t126 * t145 + t125 * qJD(3) + (-r_i_i_C(3) * t124 + (-r_i_i_C(1) * t120 - t113) * t125) * t129 + t152 * r_i_i_C(2);
	t131 = t113 * t145 + r_i_i_C(2) * t136 + t125 * t142 + (-r_i_i_C(2) * t119 * t129 - qJD(3)) * t124 + (-r_i_i_C(3) + t126) * t144 + t151 * r_i_i_C(1);
	t1 = [0, 0, 0, t135 - t142, t135; sin(qJ(1)) * t148 + t131, t131, -t145, (t121 * t144 + t124 * t143) * pkin(4) + t133, t133; -cos(qJ(1)) * t148 + t132, t132, t144, (t121 * t145 - t125 * t143) * pkin(4) + t134, t134;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end