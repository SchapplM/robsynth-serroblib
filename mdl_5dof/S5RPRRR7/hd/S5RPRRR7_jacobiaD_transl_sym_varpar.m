% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPRRR7
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
% Datum: 2019-10-24 10:46
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:38
	% EndTime: 2019-10-24 10:46:38
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:38
	% EndTime: 2019-10-24 10:46:38
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
	% StartTime: 2019-10-24 10:46:38
	% EndTime: 2019-10-24 10:46:38
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->8), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t37 = -r_i_i_C(3) - qJ(2);
	t36 = r_i_i_C(1) * cos(pkin(9)) - r_i_i_C(2) * sin(pkin(9)) + pkin(1);
	t35 = cos(qJ(1));
	t34 = sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t34 * qJD(2) + (t36 * t34 + t37 * t35) * qJD(1), -qJD(1) * t34, 0, 0, 0; t35 * qJD(2) + (t37 * t34 - t36 * t35) * qJD(1), qJD(1) * t35, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:38
	% EndTime: 2019-10-24 10:46:38
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (42->23), mult. (124->40), div. (0->0), fcn. (104->6), ass. (0->20)
	t114 = sin(pkin(9));
	t115 = cos(pkin(9));
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
	% StartTime: 2019-10-24 10:46:38
	% EndTime: 2019-10-24 10:46:38
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (156->40), mult. (227->62), div. (0->0), fcn. (186->8), ass. (0->34)
	t157 = sin(pkin(9));
	t158 = cos(pkin(9));
	t161 = cos(qJ(3));
	t183 = (r_i_i_C(3) + pkin(7) + pkin(6)) * t157 + (t161 * pkin(3) + pkin(2)) * t158 + pkin(1);
	t181 = pkin(3) * qJD(3);
	t160 = sin(qJ(1));
	t180 = t158 * t160;
	t162 = cos(qJ(1));
	t179 = t158 * t162;
	t159 = sin(qJ(3));
	t178 = t159 * t160;
	t177 = t159 * t162;
	t176 = t160 * t161;
	t175 = t161 * t162;
	t155 = qJD(3) + qJD(4);
	t156 = qJ(3) + qJ(4);
	t153 = sin(t156);
	t154 = cos(t156);
	t166 = t153 * t180 + t154 * t162;
	t169 = t153 * t160 + t154 * t179;
	t147 = t166 * qJD(1) - t169 * t155;
	t167 = t153 * t179 - t154 * t160;
	t168 = -t153 * t162 + t154 * t180;
	t148 = t168 * qJD(1) + t167 * t155;
	t174 = t147 * r_i_i_C(1) + t148 * r_i_i_C(2);
	t149 = t167 * qJD(1) + t168 * t155;
	t150 = t169 * qJD(1) - t166 * t155;
	t173 = t149 * r_i_i_C(1) + t150 * r_i_i_C(2);
	t172 = r_i_i_C(1) * t154 * t155;
	t170 = -pkin(3) * t159 - qJ(2);
	t165 = t158 * t177 - t176;
	t164 = t158 * t178 + t175;
	t151 = t157 * t155 * t153 * r_i_i_C(2);
	t1 = [0, 0, t151 + (-t161 * t181 - t172) * t157, -t157 * t172 + t151, 0; t148 * r_i_i_C(1) - t147 * r_i_i_C(2) - t160 * qJD(2) + t165 * t181 + (t183 * t160 + t170 * t162) * qJD(1), -qJD(1) * t160, ((t158 * t176 - t177) * qJD(3) + t165 * qJD(1)) * pkin(3) + t173, t173, 0; -t150 * r_i_i_C(1) + t149 * r_i_i_C(2) + t162 * qJD(2) + t164 * t181 + (t170 * t160 - t183 * t162) * qJD(1), qJD(1) * t162, ((-t158 * t175 - t178) * qJD(3) + t164 * qJD(1)) * pkin(3) + t174, t174, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:46:38
	% EndTime: 2019-10-24 10:46:39
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (375->58), mult. (334->79), div. (0->0), fcn. (269->10), ass. (0->41)
	t174 = qJ(3) + qJ(4);
	t169 = cos(t174);
	t179 = cos(qJ(3));
	t163 = t179 * pkin(3) + pkin(4) * t169;
	t175 = sin(pkin(9));
	t176 = cos(pkin(9));
	t199 = (r_i_i_C(3) + pkin(8) + pkin(7) + pkin(6)) * t175 + (pkin(2) + t163) * t176 + pkin(1);
	t172 = qJD(3) + qJD(4);
	t198 = pkin(4) * t172;
	t196 = pkin(3) * qJD(3);
	t178 = sin(qJ(1));
	t195 = t176 * t178;
	t180 = cos(qJ(1));
	t194 = t176 * t180;
	t168 = sin(t174);
	t177 = sin(qJ(3));
	t158 = -t168 * t198 - t177 * t196;
	t193 = t178 * t158;
	t192 = t180 * t158;
	t162 = t177 * pkin(3) + pkin(4) * t168;
	t191 = -qJ(2) - t162;
	t167 = qJD(5) + t172;
	t170 = qJ(5) + t174;
	t165 = sin(t170);
	t166 = cos(t170);
	t181 = t165 * t195 + t166 * t180;
	t184 = t165 * t178 + t166 * t194;
	t154 = t181 * qJD(1) - t184 * t167;
	t182 = t165 * t194 - t166 * t178;
	t183 = -t165 * t180 + t166 * t195;
	t155 = t183 * qJD(1) + t182 * t167;
	t190 = t154 * r_i_i_C(1) + t155 * r_i_i_C(2);
	t156 = t182 * qJD(1) + t183 * t167;
	t157 = t184 * qJD(1) - t181 * t167;
	t189 = t156 * r_i_i_C(1) + t157 * r_i_i_C(2);
	t187 = t169 * t198;
	t159 = t179 * t196 + t187;
	t188 = qJD(2) + t159;
	t186 = r_i_i_C(1) * t166 * t167;
	t160 = t175 * t167 * t165 * r_i_i_C(2);
	t1 = [0, 0, t160 + (-t159 - t186) * t175, t160 + (-t186 - t187) * t175, -t175 * t186 + t160; -t176 * t192 + t155 * r_i_i_C(1) - t154 * r_i_i_C(2) - t188 * t178 + (t199 * t178 + t191 * t180) * qJD(1), -qJD(1) * t178, t159 * t195 + t192 + (t162 * t194 - t163 * t178) * qJD(1) + t189, ((-t168 * t180 + t169 * t195) * t172 + (t168 * t194 - t169 * t178) * qJD(1)) * pkin(4) + t189, t189; -t176 * t193 - t157 * r_i_i_C(1) + t156 * r_i_i_C(2) + t188 * t180 + (t191 * t178 - t199 * t180) * qJD(1), qJD(1) * t180, -t159 * t194 + t193 + (t162 * t195 + t163 * t180) * qJD(1) + t190, ((-t168 * t178 - t169 * t194) * t172 + (t168 * t195 + t169 * t180) * qJD(1)) * pkin(4) + t190, t190;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end