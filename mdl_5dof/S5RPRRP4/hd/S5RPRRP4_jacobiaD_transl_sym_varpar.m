% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRP4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRP4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_jacobiaD_transl_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
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
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
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
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
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
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (156->40), mult. (227->62), div. (0->0), fcn. (186->8), ass. (0->34)
	t155 = sin(pkin(8));
	t156 = cos(pkin(8));
	t159 = cos(qJ(3));
	t181 = (r_i_i_C(3) + pkin(7) + pkin(6)) * t155 + (t159 * pkin(3) + pkin(2)) * t156 + pkin(1);
	t179 = pkin(3) * qJD(3);
	t158 = sin(qJ(1));
	t178 = t156 * t158;
	t160 = cos(qJ(1));
	t177 = t156 * t160;
	t157 = sin(qJ(3));
	t176 = t157 * t158;
	t175 = t157 * t160;
	t174 = t158 * t159;
	t173 = t159 * t160;
	t153 = qJD(3) + qJD(4);
	t154 = qJ(3) + qJ(4);
	t151 = sin(t154);
	t152 = cos(t154);
	t164 = t151 * t178 + t152 * t160;
	t167 = t151 * t158 + t152 * t177;
	t145 = t164 * qJD(1) - t167 * t153;
	t165 = -t151 * t177 + t152 * t158;
	t166 = t151 * t160 - t152 * t178;
	t146 = t166 * qJD(1) + t165 * t153;
	t172 = -t145 * r_i_i_C(1) + t146 * r_i_i_C(2);
	t147 = t165 * qJD(1) + t166 * t153;
	t148 = t167 * qJD(1) - t164 * t153;
	t171 = t147 * r_i_i_C(1) - t148 * r_i_i_C(2);
	t170 = r_i_i_C(1) * t152 * t153;
	t168 = pkin(3) * t157 + qJ(2);
	t163 = -t156 * t175 + t174;
	t162 = -t156 * t176 - t173;
	t149 = t155 * t153 * t151 * r_i_i_C(2);
	t1 = [0, 0, t149 + (-t159 * t179 - t170) * t155, -t155 * t170 + t149, 0; t146 * r_i_i_C(1) + t145 * r_i_i_C(2) + t158 * qJD(2) + t163 * t179 + (-t181 * t158 + t168 * t160) * qJD(1), qJD(1) * t158, ((-t156 * t174 + t175) * qJD(3) + t163 * qJD(1)) * pkin(3) + t171, t171, 0; t148 * r_i_i_C(1) + t147 * r_i_i_C(2) - t160 * qJD(2) + t162 * t179 + (t168 * t158 + t181 * t160) * qJD(1), -qJD(1) * t160, ((t156 * t173 + t176) * qJD(3) + t162 * qJD(1)) * pkin(3) + t172, t172, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 11:51:05
	% EndTime: 2020-01-03 11:51:05
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (226->50), mult. (292->72), div. (0->0), fcn. (233->8), ass. (0->38)
	t163 = qJ(3) + qJ(4);
	t159 = cos(t163);
	t168 = cos(qJ(3));
	t156 = t168 * pkin(3) + pkin(4) * t159;
	t164 = sin(pkin(8));
	t165 = cos(pkin(8));
	t189 = (r_i_i_C(3) + qJ(5) + pkin(7) + pkin(6)) * t164 + (pkin(2) + t156) * t165 + pkin(1);
	t162 = qJD(3) + qJD(4);
	t158 = sin(t163);
	t169 = cos(qJ(1));
	t167 = sin(qJ(1));
	t183 = t165 * t167;
	t170 = t158 * t183 + t159 * t169;
	t182 = t165 * t169;
	t173 = t158 * t167 + t159 * t182;
	t147 = t170 * qJD(1) - t173 * t162;
	t188 = pkin(4) * t158;
	t186 = pkin(3) * qJD(3);
	t185 = t159 * t162;
	t184 = t162 * t164;
	t166 = sin(qJ(3));
	t155 = t166 * pkin(3) + t188;
	t181 = qJ(2) + t155;
	t171 = -t158 * t182 + t159 * t167;
	t172 = t158 * t169 - t159 * t183;
	t148 = t172 * qJD(1) + t171 * t162;
	t180 = -t147 * r_i_i_C(1) + t148 * r_i_i_C(2);
	t149 = t171 * qJD(1) + t172 * t162;
	t150 = t173 * qJD(1) - t170 * t162;
	t179 = t149 * r_i_i_C(1) - t150 * r_i_i_C(2);
	t178 = qJD(1) * t167;
	t177 = qJD(1) * t169;
	t152 = pkin(4) * t185 + t168 * t186;
	t176 = qJD(2) + t152;
	t151 = -t162 * t188 - t166 * t186;
	t174 = qJD(5) * t164 + t151 * t165;
	t153 = t158 * r_i_i_C(2) * t184;
	t1 = [0, 0, t153 + (-r_i_i_C(1) * t185 - t152) * t164, t153 + (-pkin(4) - r_i_i_C(1)) * t159 * t184, 0; t148 * r_i_i_C(1) + t147 * r_i_i_C(2) + t174 * t169 + t176 * t167 + (-t189 * t167 + t181 * t169) * qJD(1), t178, -t152 * t183 - t169 * t151 + (-t155 * t182 + t156 * t167) * qJD(1) + t179, t149 * pkin(4) + t179, t164 * t177; t150 * r_i_i_C(1) + t149 * r_i_i_C(2) - t176 * t169 + t174 * t167 + (t181 * t167 + t189 * t169) * qJD(1), -t177, t152 * t182 - t167 * t151 + (-t155 * t183 - t156 * t169) * qJD(1) + t180, -t147 * pkin(4) + t180, t164 * t178;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end