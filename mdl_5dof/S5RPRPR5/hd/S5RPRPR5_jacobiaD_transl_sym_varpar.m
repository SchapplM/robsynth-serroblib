% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:44:05
	% EndTime: 2020-01-03 11:44:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:44:05
	% EndTime: 2020-01-03 11:44:05
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
	% StartTime: 2020-01-03 11:44:05
	% EndTime: 2020-01-03 11:44:05
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->8), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t37 = r_i_i_C(3) + qJ(2);
	t36 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(1);
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t34 * qJD(2) + (-t36 * t34 + t37 * t35) * qJD(1), qJD(1) * t34, 0, 0, 0; -t35 * qJD(2) + (t37 * t34 + t36 * t35) * qJD(1), -qJD(1) * t35, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:44:05
	% EndTime: 2020-01-03 11:44:05
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (42->23), mult. (124->40), div. (0->0), fcn. (104->6), ass. (0->20)
	t110 = sin(pkin(8));
	t111 = cos(pkin(8));
	t126 = pkin(2) * t111 + pkin(1) + (pkin(6) + r_i_i_C(3)) * t110;
	t112 = sin(qJ(3));
	t113 = sin(qJ(1));
	t124 = t112 * t113;
	t115 = cos(qJ(1));
	t123 = t112 * t115;
	t114 = cos(qJ(3));
	t122 = t113 * t114;
	t121 = t114 * t115;
	t119 = t111 * t121 + t124;
	t118 = -t111 * t122 + t123;
	t117 = -t111 * t123 + t122;
	t116 = t111 * t124 + t121;
	t109 = t119 * qJD(1) - t116 * qJD(3);
	t108 = t117 * qJD(1) + t118 * qJD(3);
	t107 = t118 * qJD(1) + t117 * qJD(3);
	t106 = t116 * qJD(1) - t119 * qJD(3);
	t1 = [0, 0, (-r_i_i_C(1) * t114 + r_i_i_C(2) * t112) * t110 * qJD(3), 0, 0; t107 * r_i_i_C(1) + t106 * r_i_i_C(2) + t113 * qJD(2) + (qJ(2) * t115 - t126 * t113) * qJD(1), qJD(1) * t113, t108 * r_i_i_C(1) - t109 * r_i_i_C(2), 0, 0; t109 * r_i_i_C(1) + t108 * r_i_i_C(2) - t115 * qJD(2) + (qJ(2) * t113 + t126 * t115) * qJD(1), -qJD(1) * t115, -t106 * r_i_i_C(1) + t107 * r_i_i_C(2), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:44:05
	% EndTime: 2020-01-03 11:44:05
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (97->40), mult. (185->66), div. (0->0), fcn. (150->8), ass. (0->33)
	t122 = sin(pkin(8));
	t123 = cos(pkin(8));
	t127 = cos(qJ(3));
	t148 = t127 * pkin(3);
	t149 = (r_i_i_C(3) + qJ(4) + pkin(6)) * t122 + (pkin(2) + t148) * t123 + pkin(1);
	t146 = pkin(3) * qJD(3);
	t126 = sin(qJ(1));
	t145 = t123 * t126;
	t128 = cos(qJ(1));
	t144 = t123 * t128;
	t125 = sin(qJ(3));
	t143 = t125 * t126;
	t142 = t125 * t128;
	t141 = t126 * t127;
	t140 = t127 * t128;
	t139 = qJD(1) * t126;
	t138 = qJD(1) * t128;
	t137 = qJD(4) * t122;
	t135 = pkin(3) * t125 + qJ(2);
	t121 = qJ(3) + pkin(9);
	t119 = sin(t121);
	t120 = cos(t121);
	t134 = t119 * t126 + t120 * t144;
	t133 = t119 * t128 - t120 * t145;
	t132 = -t119 * t144 + t120 * t126;
	t131 = t119 * t145 + t120 * t128;
	t130 = -t123 * t142 + t141;
	t129 = -t123 * t143 - t140;
	t117 = t134 * qJD(1) - t131 * qJD(3);
	t116 = t132 * qJD(1) + t133 * qJD(3);
	t115 = t133 * qJD(1) + t132 * qJD(3);
	t114 = t131 * qJD(1) - t134 * qJD(3);
	t1 = [0, 0, (-r_i_i_C(1) * t120 + r_i_i_C(2) * t119 - t148) * t122 * qJD(3), 0, 0; t128 * t137 + t115 * r_i_i_C(1) + t114 * r_i_i_C(2) + t126 * qJD(2) + t130 * t146 + (-t149 * t126 + t135 * t128) * qJD(1), t139, t116 * r_i_i_C(1) - t117 * r_i_i_C(2) + ((-t123 * t141 + t142) * qJD(3) + t130 * qJD(1)) * pkin(3), t122 * t138, 0; t126 * t137 + t117 * r_i_i_C(1) + t116 * r_i_i_C(2) - t128 * qJD(2) + t129 * t146 + (t135 * t126 + t149 * t128) * qJD(1), -t138, -t114 * r_i_i_C(1) + t115 * r_i_i_C(2) + ((t123 * t140 + t143) * qJD(3) + t129 * qJD(1)) * pkin(3), t122 * t139, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:44:05
	% EndTime: 2020-01-03 11:44:05
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (246->46), mult. (261->64), div. (0->0), fcn. (211->10), ass. (0->34)
	t167 = qJ(3) + pkin(9);
	t157 = pkin(4) * cos(t167) + cos(qJ(3)) * pkin(3);
	t168 = sin(pkin(8));
	t169 = cos(pkin(8));
	t190 = (r_i_i_C(3) + pkin(7) + qJ(4) + pkin(6)) * t168 + (pkin(2) + t157) * t169 + pkin(1);
	t171 = sin(qJ(1));
	t188 = t169 * t171;
	t173 = cos(qJ(1));
	t187 = t169 * t173;
	t156 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t167);
	t186 = qJ(2) + t156;
	t166 = qJD(3) + qJD(5);
	t163 = qJ(5) + t167;
	t159 = sin(t163);
	t160 = cos(t163);
	t174 = t159 * t188 + t160 * t173;
	t177 = t159 * t171 + t160 * t187;
	t148 = t174 * qJD(1) - t177 * t166;
	t175 = -t159 * t187 + t160 * t171;
	t176 = t159 * t173 - t160 * t188;
	t149 = t176 * qJD(1) + t175 * t166;
	t185 = -t148 * r_i_i_C(1) + t149 * r_i_i_C(2);
	t150 = t175 * qJD(1) + t176 * t166;
	t151 = t177 * qJD(1) - t174 * t166;
	t184 = t150 * r_i_i_C(1) - t151 * r_i_i_C(2);
	t183 = qJD(1) * t171;
	t182 = qJD(1) * t173;
	t153 = t157 * qJD(3);
	t181 = qJD(2) + t153;
	t180 = r_i_i_C(1) * t160 * t166;
	t152 = t156 * qJD(3);
	t178 = qJD(4) * t168 - t152 * t169;
	t154 = t168 * t166 * t159 * r_i_i_C(2);
	t1 = [0, 0, t154 + (-t153 - t180) * t168, 0, -t168 * t180 + t154; t149 * r_i_i_C(1) + t148 * r_i_i_C(2) + t178 * t173 + t181 * t171 + (-t190 * t171 + t186 * t173) * qJD(1), t183, -t153 * t188 + t173 * t152 + (-t156 * t187 + t157 * t171) * qJD(1) + t184, t168 * t182, t184; t151 * r_i_i_C(1) + t150 * r_i_i_C(2) - t181 * t173 + t178 * t171 + (t186 * t171 + t190 * t173) * qJD(1), -t182, t153 * t187 + t171 * t152 + (-t156 * t188 - t157 * t173) * qJD(1) + t185, t168 * t183, t185;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end