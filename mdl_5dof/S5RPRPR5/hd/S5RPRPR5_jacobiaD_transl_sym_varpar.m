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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
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
	% StartTime: 2022-01-23 09:26:46
	% EndTime: 2022-01-23 09:26:46
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
	% StartTime: 2022-01-23 09:26:47
	% EndTime: 2022-01-23 09:26:47
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (41->22), mult. (116->40), div. (0->0), fcn. (100->6), ass. (0->20)
	t143 = sin(qJ(3));
	t144 = sin(qJ(1));
	t155 = t143 * t144;
	t146 = cos(qJ(1));
	t154 = t143 * t146;
	t145 = cos(qJ(3));
	t153 = t144 * t145;
	t152 = t145 * t146;
	t141 = sin(pkin(8));
	t142 = cos(pkin(8));
	t151 = -pkin(2) * t142 - pkin(1) + (-pkin(6) - r_i_i_C(3)) * t141;
	t150 = -t142 * t152 - t155;
	t149 = t142 * t153 - t154;
	t148 = t142 * t154 - t153;
	t147 = t142 * t155 + t152;
	t139 = qJD(1) * t150 + qJD(3) * t147;
	t138 = qJD(1) * t148 + qJD(3) * t149;
	t137 = qJD(1) * t149 + qJD(3) * t148;
	t136 = qJD(1) * t147 + qJD(3) * t150;
	t1 = [t139 * r_i_i_C(1) + t138 * r_i_i_C(2) + t146 * qJD(2) + (-qJ(2) * t144 + t146 * t151) * qJD(1), qJD(1) * t146, r_i_i_C(1) * t136 + r_i_i_C(2) * t137, 0, 0; -t137 * r_i_i_C(1) + t136 * r_i_i_C(2) + t144 * qJD(2) + (qJ(2) * t146 + t144 * t151) * qJD(1), qJD(1) * t144, -r_i_i_C(1) * t138 + r_i_i_C(2) * t139, 0, 0; 0, 0, (-r_i_i_C(1) * t145 + r_i_i_C(2) * t143) * t141 * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:47
	% EndTime: 2022-01-23 09:26:47
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (97->40), mult. (171->67), div. (0->0), fcn. (142->8), ass. (0->32)
	t157 = cos(qJ(3));
	t173 = t157 * pkin(3);
	t172 = pkin(3) * qJD(3);
	t154 = cos(pkin(8));
	t156 = sin(qJ(1));
	t171 = t154 * t156;
	t158 = cos(qJ(1));
	t170 = t154 * t158;
	t155 = sin(qJ(3));
	t169 = t155 * t156;
	t168 = t155 * t158;
	t167 = t156 * t157;
	t166 = t157 * t158;
	t165 = qJD(1) * t156;
	t164 = qJD(1) * t158;
	t153 = sin(pkin(8));
	t163 = -pkin(1) - (pkin(2) + t173) * t154 + (-r_i_i_C(3) - qJ(4) - pkin(6)) * t153;
	t152 = qJ(3) + pkin(9);
	t150 = sin(t152);
	t151 = cos(t152);
	t162 = -t150 * t156 - t151 * t170;
	t161 = -t150 * t158 + t151 * t171;
	t160 = t150 * t170 - t151 * t156;
	t159 = t150 * t171 + t151 * t158;
	t149 = t155 * pkin(3) + qJ(2);
	t148 = t157 * t172 + qJD(2);
	t147 = -t154 * t155 * t172 + qJD(4) * t153;
	t145 = t162 * qJD(1) + t159 * qJD(3);
	t144 = t160 * qJD(1) + t161 * qJD(3);
	t143 = t161 * qJD(1) + t160 * qJD(3);
	t142 = t159 * qJD(1) + t162 * qJD(3);
	t1 = [t145 * r_i_i_C(1) + t144 * r_i_i_C(2) - t147 * t156 + t148 * t158 + (-t149 * t156 + t163 * t158) * qJD(1), t164, t142 * r_i_i_C(1) + t143 * r_i_i_C(2) + ((-t154 * t166 - t169) * qJD(3) + (t154 * t169 + t166) * qJD(1)) * pkin(3), -t153 * t165, 0; -t143 * r_i_i_C(1) + t142 * r_i_i_C(2) + t147 * t158 + t148 * t156 + (t149 * t158 + t163 * t156) * qJD(1), t165, -t144 * r_i_i_C(1) + t145 * r_i_i_C(2) + ((-t154 * t167 + t168) * qJD(3) + (-t154 * t168 + t167) * qJD(1)) * pkin(3), t153 * t164, 0; 0, 0, (-r_i_i_C(1) * t151 + r_i_i_C(2) * t150 - t173) * t153 * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:26:47
	% EndTime: 2022-01-23 09:26:47
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (246->46), mult. (261->64), div. (0->0), fcn. (211->10), ass. (0->34)
	t187 = qJ(3) + pkin(9);
	t177 = pkin(4) * cos(t187) + cos(qJ(3)) * pkin(3);
	t189 = cos(pkin(8));
	t191 = sin(qJ(1));
	t208 = t189 * t191;
	t193 = cos(qJ(1));
	t207 = t189 * t193;
	t176 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t187);
	t206 = qJ(2) + t176;
	t186 = qJD(3) + qJD(5);
	t183 = qJ(5) + t187;
	t179 = sin(t183);
	t180 = cos(t183);
	t195 = t179 * t208 + t180 * t193;
	t198 = -t179 * t191 - t180 * t207;
	t168 = t195 * qJD(1) + t198 * t186;
	t196 = t179 * t207 - t180 * t191;
	t197 = -t179 * t193 + t180 * t208;
	t169 = t197 * qJD(1) + t196 * t186;
	t205 = t168 * r_i_i_C(1) + t169 * r_i_i_C(2);
	t170 = t196 * qJD(1) + t197 * t186;
	t171 = t198 * qJD(1) + t195 * t186;
	t204 = -t170 * r_i_i_C(1) + t171 * r_i_i_C(2);
	t203 = qJD(1) * t191;
	t202 = qJD(1) * t193;
	t173 = t177 * qJD(3);
	t201 = qJD(2) + t173;
	t200 = r_i_i_C(1) * t180 * t186;
	t172 = t176 * qJD(3);
	t188 = sin(pkin(8));
	t199 = qJD(4) * t188 - t172 * t189;
	t194 = -(pkin(2) + t177) * t189 - pkin(1) + (-r_i_i_C(3) - pkin(7) - qJ(4) - pkin(6)) * t188;
	t174 = t188 * t186 * t179 * r_i_i_C(2);
	t1 = [t171 * r_i_i_C(1) + t170 * r_i_i_C(2) + t201 * t193 - t199 * t191 + (-t206 * t191 + t194 * t193) * qJD(1), t202, -t173 * t207 - t191 * t172 + (t176 * t208 + t177 * t193) * qJD(1) + t205, -t188 * t203, t205; -t169 * r_i_i_C(1) + t168 * r_i_i_C(2) + t199 * t193 + t201 * t191 + (t194 * t191 + t206 * t193) * qJD(1), t203, -t173 * t208 + t193 * t172 + (-t176 * t207 + t177 * t191) * qJD(1) + t204, t188 * t202, t204; 0, 0, t174 + (-t173 - t200) * t188, 0, -t188 * t200 + t174;];
	JaD_transl = t1;
end