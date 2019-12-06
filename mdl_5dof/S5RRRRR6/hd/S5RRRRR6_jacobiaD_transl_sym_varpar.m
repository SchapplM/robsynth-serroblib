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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
	% StartTime: 2019-12-05 19:01:29
	% EndTime: 2019-12-05 19:01:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 19:01:28
	% EndTime: 2019-12-05 19:01:28
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
	% StartTime: 2019-12-05 19:01:29
	% EndTime: 2019-12-05 19:01:29
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
	% StartTime: 2019-12-05 19:01:29
	% EndTime: 2019-12-05 19:01:29
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (93->20), mult. (104->34), div. (0->0), fcn. (64->6), ass. (0->19)
	t89 = sin(qJ(3));
	t100 = qJD(3) * t89;
	t87 = qJD(1) + qJD(2);
	t90 = cos(qJ(3));
	t102 = t87 * t90;
	t88 = qJ(1) + qJ(2);
	t85 = sin(t88);
	t86 = cos(t88);
	t106 = t86 * t100 + t85 * t102;
	t103 = t87 * t89;
	t99 = qJD(3) * t90;
	t105 = t86 * t103 + t85 * t99;
	t104 = -pkin(7) - r_i_i_C(3);
	t101 = pkin(1) * qJD(1);
	t96 = t85 * t100;
	t93 = t86 * t99;
	t92 = r_i_i_C(2) * t93 + (t104 * t86 + (-r_i_i_C(2) * t89 + pkin(2)) * t85) * t87 + t106 * r_i_i_C(1);
	t91 = r_i_i_C(1) * t96 + ((-r_i_i_C(1) * t90 - pkin(2)) * t86 + t104 * t85) * t87 + t105 * r_i_i_C(2);
	t1 = [0, 0, (-r_i_i_C(1) * t89 - r_i_i_C(2) * t90) * qJD(3), 0, 0; sin(qJ(1)) * t101 + t92, t92, (t86 * t102 - t96) * r_i_i_C(2) + t105 * r_i_i_C(1), 0, 0; -cos(qJ(1)) * t101 + t91, t91, t106 * r_i_i_C(2) + (t85 * t103 - t93) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 19:01:29
	% EndTime: 2019-12-05 19:01:29
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (195->30), mult. (162->43), div. (0->0), fcn. (103->8), ass. (0->32)
	t124 = qJ(1) + qJ(2);
	t118 = sin(t124);
	t122 = qJD(1) + qJD(2);
	t145 = t118 * t122;
	t120 = cos(t124);
	t143 = t120 * t122;
	t123 = qJ(3) + qJ(4);
	t117 = sin(t123);
	t119 = cos(t123);
	t121 = qJD(3) + qJD(4);
	t146 = t118 * t121;
	t152 = t117 * t143 + t119 * t146;
	t144 = t120 * t121;
	t151 = t117 * t144 + t119 * t145;
	t125 = sin(qJ(3));
	t139 = qJD(3) * t125 * pkin(3);
	t150 = t139 + (-r_i_i_C(3) - pkin(8) - pkin(7)) * t122;
	t149 = r_i_i_C(1) * t117;
	t148 = r_i_i_C(2) * t119;
	t147 = pkin(1) * qJD(1);
	t142 = t122 * t125;
	t126 = cos(qJ(3));
	t140 = qJD(3) * t126;
	t138 = t117 * t146;
	t133 = t119 * t144;
	t132 = (-t148 - t149) * t121;
	t131 = -r_i_i_C(1) * t133 + r_i_i_C(2) * t151 + t145 * t149;
	t130 = r_i_i_C(1) * t152 - r_i_i_C(2) * t138 + t143 * t148;
	t116 = t126 * pkin(3) + pkin(2);
	t129 = t116 * t145 + (-t117 * t145 + t133) * r_i_i_C(2) + t150 * t120 + t151 * r_i_i_C(1);
	t128 = r_i_i_C(1) * t138 + (-r_i_i_C(1) * t119 - t116) * t143 + t150 * t118 + t152 * r_i_i_C(2);
	t1 = [0, 0, t132 - t139, t132, 0; sin(qJ(1)) * t147 + t129, t129, (t118 * t140 + t120 * t142) * pkin(3) + t130, t130, 0; -cos(qJ(1)) * t147 + t128, t128, (t118 * t142 - t120 * t140) * pkin(3) + t131, t131, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 19:01:29
	% EndTime: 2019-12-05 19:01:29
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (358->44), mult. (224->56), div. (0->0), fcn. (143->10), ass. (0->44)
	t134 = qJ(1) + qJ(2);
	t128 = cos(t134);
	t131 = qJD(1) + qJD(2);
	t155 = t128 * t131;
	t133 = qJ(3) + qJ(4);
	t129 = qJ(5) + t133;
	t122 = sin(t129);
	t123 = cos(t129);
	t130 = qJD(3) + qJD(4);
	t124 = qJD(5) + t130;
	t126 = sin(t134);
	t159 = t124 * t126;
	t167 = t122 * t155 + t123 * t159;
	t157 = t126 * t131;
	t158 = t124 * t128;
	t166 = t122 * t158 + t123 * t157;
	t135 = sin(qJ(3));
	t160 = pkin(3) * qJD(3);
	t152 = t135 * t160;
	t125 = sin(t133);
	t164 = pkin(4) * t125;
	t153 = t130 * t164;
	t165 = t152 + t153 + (-r_i_i_C(3) - pkin(9) - pkin(8) - pkin(7)) * t131;
	t163 = r_i_i_C(1) * t123;
	t162 = r_i_i_C(2) * t123;
	t161 = pkin(1) * qJD(1);
	t127 = cos(t133);
	t156 = t127 * t130;
	t151 = t122 * t159;
	t149 = t122 * t157;
	t146 = t123 * t158;
	t144 = t167 * r_i_i_C(1) + t155 * t162;
	t143 = r_i_i_C(1) * t149 + t166 * r_i_i_C(2);
	t142 = (-r_i_i_C(1) * t122 - t162) * t124;
	t141 = -r_i_i_C(1) * t146 + t143;
	t140 = -r_i_i_C(2) * t151 + t144;
	t139 = t142 - t153;
	t136 = cos(qJ(3));
	t120 = t136 * pkin(3) + pkin(4) * t127 + pkin(2);
	t138 = t120 * t157 + t165 * t128 + (-t149 + t146) * r_i_i_C(2) + t166 * r_i_i_C(1);
	t137 = r_i_i_C(1) * t151 + (-t120 - t163) * t155 + t165 * t126 + t167 * r_i_i_C(2);
	t121 = -t135 * pkin(3) - t164;
	t111 = -pkin(4) * t156 - t136 * t160;
	t1 = [0, 0, t139 - t152, t139, t142; sin(qJ(1)) * t161 + t138, t138, -t121 * t155 + (-r_i_i_C(2) * t122 * t124 - t111) * t126 + t144, (t125 * t155 + t126 * t156) * pkin(4) + t140, t140; -cos(qJ(1)) * t161 + t137, t137, -t121 * t157 + (-t124 * t163 + t111) * t128 + t143, (t125 * t157 - t128 * t156) * pkin(4) + t141, t141;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end