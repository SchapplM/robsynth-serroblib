% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPPRR4_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPPRR4_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPRR4_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:41
	% EndTime: 2022-01-23 09:17:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->7), mult. (28->12), div. (0->0), fcn. (18->4), ass. (0->5)
	t16 = r_i_i_C(3) + qJ(2);
	t15 = -r_i_i_C(1) * cos(pkin(8)) + r_i_i_C(2) * sin(pkin(8)) - pkin(1);
	t14 = cos(qJ(1));
	t13 = sin(qJ(1));
	t1 = [t14 * qJD(2) + (-t16 * t13 + t15 * t14) * qJD(1), qJD(1) * t14, 0, 0, 0; t13 * qJD(2) + (t15 * t13 + t16 * t14) * qJD(1), qJD(1) * t13, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:42
	% EndTime: 2022-01-23 09:17:42
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (21->14), mult. (54->21), div. (0->0), fcn. (44->6), ass. (0->11)
	t118 = sin(qJ(1));
	t124 = qJD(1) * t118;
	t119 = cos(qJ(1));
	t123 = qJD(1) * t119;
	t115 = sin(pkin(8));
	t122 = t115 * qJD(3);
	t114 = sin(pkin(9));
	t116 = cos(pkin(9));
	t121 = r_i_i_C(1) * t114 + r_i_i_C(2) * t116 + qJ(2);
	t120 = -pkin(1) + (-r_i_i_C(1) * t116 + r_i_i_C(2) * t114 - pkin(2)) * cos(pkin(8)) + (-r_i_i_C(3) - qJ(3)) * t115;
	t1 = [-t118 * t122 + t119 * qJD(2) + (-t118 * t121 + t119 * t120) * qJD(1), t123, -t115 * t124, 0, 0; t119 * t122 + t118 * qJD(2) + (t118 * t120 + t119 * t121) * qJD(1), t124, t115 * t123, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:42
	% EndTime: 2022-01-23 09:17:42
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (84->29), mult. (128->49), div. (0->0), fcn. (112->8), ass. (0->23)
	t148 = cos(pkin(8));
	t149 = sin(qJ(1));
	t160 = t148 * t149;
	t150 = cos(qJ(1));
	t159 = t148 * t150;
	t158 = qJD(1) * t149;
	t157 = qJD(1) * t150;
	t147 = sin(pkin(8));
	t156 = qJD(3) * t147;
	t155 = -(cos(pkin(9)) * pkin(3) + pkin(2)) * t148 - pkin(1) + (-r_i_i_C(3) - qJ(3) - pkin(6)) * t147;
	t146 = pkin(9) + qJ(4);
	t144 = sin(t146);
	t145 = cos(t146);
	t154 = -t144 * t149 - t145 * t159;
	t153 = -t144 * t150 + t145 * t160;
	t152 = t144 * t159 - t145 * t149;
	t151 = t144 * t160 + t145 * t150;
	t143 = sin(pkin(9)) * pkin(3) + qJ(2);
	t141 = t154 * qJD(1) + t151 * qJD(4);
	t140 = t152 * qJD(1) + t153 * qJD(4);
	t139 = t153 * qJD(1) + t152 * qJD(4);
	t138 = t151 * qJD(1) + t154 * qJD(4);
	t1 = [-t149 * t156 + t141 * r_i_i_C(1) + t140 * r_i_i_C(2) + qJD(2) * t150 + (-t143 * t149 + t155 * t150) * qJD(1), t157, -t147 * t158, t138 * r_i_i_C(1) + t139 * r_i_i_C(2), 0; t150 * t156 - t139 * r_i_i_C(1) + t138 * r_i_i_C(2) + qJD(2) * t149 + (t143 * t150 + t155 * t149) * qJD(1), t158, t147 * t157, -t140 * r_i_i_C(1) + t141 * r_i_i_C(2), 0; 0, 0, 0, (-r_i_i_C(1) * t145 + r_i_i_C(2) * t144) * t147 * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:17:42
	% EndTime: 2022-01-23 09:17:42
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (233->46), mult. (239->69), div. (0->0), fcn. (198->10), ass. (0->34)
	t203 = pkin(4) * qJD(4);
	t184 = cos(pkin(8));
	t185 = sin(qJ(1));
	t202 = t184 * t185;
	t186 = cos(qJ(1));
	t201 = t184 * t186;
	t181 = pkin(9) + qJ(4);
	t177 = sin(t181);
	t200 = qJ(2) + pkin(4) * t177 + sin(pkin(9)) * pkin(3);
	t182 = qJD(4) + qJD(5);
	t179 = qJ(5) + t181;
	t175 = sin(t179);
	t176 = cos(t179);
	t190 = t175 * t202 + t176 * t186;
	t193 = -t175 * t185 - t176 * t201;
	t168 = t190 * qJD(1) + t193 * t182;
	t191 = t175 * t201 - t176 * t185;
	t192 = -t175 * t186 + t176 * t202;
	t169 = t192 * qJD(1) + t191 * t182;
	t199 = t168 * r_i_i_C(1) + t169 * r_i_i_C(2);
	t170 = t191 * qJD(1) + t192 * t182;
	t171 = t193 * qJD(1) + t190 * t182;
	t198 = -t170 * r_i_i_C(1) + t171 * r_i_i_C(2);
	t197 = qJD(1) * t185;
	t196 = qJD(1) * t186;
	t183 = sin(pkin(8));
	t195 = qJD(3) * t183;
	t194 = r_i_i_C(1) * t176 * t182;
	t178 = cos(t181);
	t189 = -t177 * t201 + t178 * t185;
	t188 = t177 * t202 + t178 * t186;
	t187 = -(pkin(4) * t178 + cos(pkin(9)) * pkin(3) + pkin(2)) * t184 - pkin(1) + (-r_i_i_C(3) - pkin(7) - pkin(6) - qJ(3)) * t183;
	t172 = t183 * t182 * t175 * r_i_i_C(2);
	t1 = [-t185 * t195 + t171 * r_i_i_C(1) + t170 * r_i_i_C(2) + t186 * qJD(2) + t188 * t203 + (-t200 * t185 + t187 * t186) * qJD(1), t196, -t183 * t197, ((-t177 * t185 - t178 * t201) * qJD(4) + t188 * qJD(1)) * pkin(4) + t199, t199; t186 * t195 - t169 * r_i_i_C(1) + t168 * r_i_i_C(2) + t185 * qJD(2) + t189 * t203 + (t187 * t185 + t200 * t186) * qJD(1), t197, t183 * t196, ((t177 * t186 - t178 * t202) * qJD(4) + t189 * qJD(1)) * pkin(4) + t198, t198; 0, 0, 0, t172 + (-t178 * t203 - t194) * t183, -t183 * t194 + t172;];
	JaD_transl = t1;
end