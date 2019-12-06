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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
	% StartTime: 2019-12-05 17:58:10
	% EndTime: 2019-12-05 17:58:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:58:10
	% EndTime: 2019-12-05 17:58:10
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
	% StartTime: 2019-12-05 17:58:10
	% EndTime: 2019-12-05 17:58:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->8), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t37 = -r_i_i_C(3) - qJ(2);
	t36 = r_i_i_C(1) * cos(pkin(8)) - r_i_i_C(2) * sin(pkin(8)) + pkin(1);
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t34 * qJD(2) + (t36 * t34 + t37 * t35) * qJD(1), -qJD(1) * t34, 0, 0, 0; t35 * qJD(2) + (t37 * t34 - t36 * t35) * qJD(1), qJD(1) * t35, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:58:10
	% EndTime: 2019-12-05 17:58:10
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (42->23), mult. (124->40), div. (0->0), fcn. (104->6), ass. (0->20)
	t114 = sin(pkin(8));
	t115 = cos(pkin(8));
	t130 = pkin(2) * t115 + pkin(1) + (pkin(6) + r_i_i_C(3)) * t114;
	t116 = sin(qJ(3));
	t117 = sin(qJ(1));
	t128 = t116 * t117;
	t119 = cos(qJ(1));
	t127 = t116 * t119;
	t118 = cos(qJ(3));
	t126 = t117 * t118;
	t125 = t118 * t119;
	t123 = t115 * t125 + t128;
	t122 = t115 * t126 - t127;
	t121 = t115 * t127 - t126;
	t120 = t115 * t128 + t125;
	t113 = t123 * qJD(1) - t120 * qJD(3);
	t112 = t121 * qJD(1) + t122 * qJD(3);
	t111 = t122 * qJD(1) + t121 * qJD(3);
	t110 = t120 * qJD(1) - t123 * qJD(3);
	t1 = [0, 0, (-r_i_i_C(1) * t118 + r_i_i_C(2) * t116) * t114 * qJD(3), 0, 0; t111 * r_i_i_C(1) - t110 * r_i_i_C(2) - t117 * qJD(2) + (-qJ(2) * t119 + t130 * t117) * qJD(1), -qJD(1) * t117, t112 * r_i_i_C(1) + t113 * r_i_i_C(2), 0, 0; -t113 * r_i_i_C(1) + t112 * r_i_i_C(2) + t119 * qJD(2) + (-qJ(2) * t117 - t130 * t119) * qJD(1), qJD(1) * t119, t110 * r_i_i_C(1) + t111 * r_i_i_C(2), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:58:10
	% EndTime: 2019-12-05 17:58:11
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (99->42), mult. (185->66), div. (0->0), fcn. (150->8), ass. (0->33)
	t126 = sin(pkin(8));
	t127 = cos(pkin(8));
	t131 = cos(qJ(3));
	t152 = t131 * pkin(3);
	t153 = (r_i_i_C(3) + qJ(4) + pkin(6)) * t126 + (pkin(2) + t152) * t127 + pkin(1);
	t150 = pkin(3) * qJD(3);
	t130 = sin(qJ(1));
	t149 = t127 * t130;
	t132 = cos(qJ(1));
	t148 = t127 * t132;
	t129 = sin(qJ(3));
	t147 = t129 * t130;
	t146 = t129 * t132;
	t145 = t130 * t131;
	t144 = t131 * t132;
	t143 = qJD(1) * t130;
	t142 = qJD(1) * t132;
	t141 = qJD(4) * t126;
	t139 = -pkin(3) * t129 - qJ(2);
	t125 = qJ(3) + pkin(9);
	t123 = sin(t125);
	t124 = cos(t125);
	t138 = t123 * t130 + t124 * t148;
	t137 = -t123 * t132 + t124 * t149;
	t136 = t123 * t148 - t124 * t130;
	t135 = t123 * t149 + t124 * t132;
	t134 = t127 * t146 - t145;
	t133 = t127 * t147 + t144;
	t121 = t138 * qJD(1) - t135 * qJD(3);
	t120 = t136 * qJD(1) + t137 * qJD(3);
	t119 = t137 * qJD(1) + t136 * qJD(3);
	t118 = t135 * qJD(1) - t138 * qJD(3);
	t1 = [0, 0, (-r_i_i_C(1) * t124 + r_i_i_C(2) * t123 - t152) * t126 * qJD(3), 0, 0; -t132 * t141 + t119 * r_i_i_C(1) - t118 * r_i_i_C(2) - t130 * qJD(2) + t134 * t150 + (t153 * t130 + t139 * t132) * qJD(1), -t143, t120 * r_i_i_C(1) + t121 * r_i_i_C(2) + ((t127 * t145 - t146) * qJD(3) + t134 * qJD(1)) * pkin(3), -t126 * t142, 0; -t130 * t141 - t121 * r_i_i_C(1) + t120 * r_i_i_C(2) + t132 * qJD(2) + t133 * t150 + (t139 * t130 - t153 * t132) * qJD(1), t142, t118 * r_i_i_C(1) + t119 * r_i_i_C(2) + ((-t127 * t144 - t147) * qJD(3) + t133 * qJD(1)) * pkin(3), -t126 * t143, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:58:10
	% EndTime: 2019-12-05 17:58:11
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (248->48), mult. (261->64), div. (0->0), fcn. (211->10), ass. (0->34)
	t169 = qJ(3) + pkin(9);
	t159 = pkin(4) * cos(t169) + cos(qJ(3)) * pkin(3);
	t170 = sin(pkin(8));
	t171 = cos(pkin(8));
	t192 = (r_i_i_C(3) + pkin(7) + qJ(4) + pkin(6)) * t170 + (pkin(2) + t159) * t171 + pkin(1);
	t173 = sin(qJ(1));
	t190 = t171 * t173;
	t175 = cos(qJ(1));
	t189 = t171 * t175;
	t158 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t169);
	t188 = -qJ(2) - t158;
	t168 = qJD(3) + qJD(5);
	t165 = qJ(5) + t169;
	t161 = sin(t165);
	t162 = cos(t165);
	t176 = t161 * t190 + t162 * t175;
	t179 = t161 * t173 + t162 * t189;
	t150 = qJD(1) * t176 - t168 * t179;
	t177 = t161 * t189 - t162 * t173;
	t178 = -t161 * t175 + t162 * t190;
	t151 = qJD(1) * t178 + t168 * t177;
	t187 = t150 * r_i_i_C(1) + t151 * r_i_i_C(2);
	t152 = qJD(1) * t177 + t168 * t178;
	t153 = qJD(1) * t179 - t168 * t176;
	t186 = t152 * r_i_i_C(1) + t153 * r_i_i_C(2);
	t185 = qJD(1) * t173;
	t184 = qJD(1) * t175;
	t155 = t159 * qJD(3);
	t183 = qJD(2) + t155;
	t182 = r_i_i_C(1) * t162 * t168;
	t154 = t158 * qJD(3);
	t180 = -qJD(4) * t170 + t154 * t171;
	t156 = t170 * t168 * t161 * r_i_i_C(2);
	t1 = [0, 0, t156 + (-t155 - t182) * t170, 0, -t170 * t182 + t156; t151 * r_i_i_C(1) - t150 * r_i_i_C(2) + t180 * t175 - t183 * t173 + (t192 * t173 + t188 * t175) * qJD(1), -t185, t155 * t190 - t175 * t154 + (t158 * t189 - t159 * t173) * qJD(1) + t186, -t170 * t184, t186; -t153 * r_i_i_C(1) + t152 * r_i_i_C(2) + t183 * t175 + t180 * t173 + (t188 * t173 - t192 * t175) * qJD(1), t184, -t155 * t189 - t173 * t154 + (t158 * t190 + t159 * t175) * qJD(1) + t187, -t170 * t185, t187;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end