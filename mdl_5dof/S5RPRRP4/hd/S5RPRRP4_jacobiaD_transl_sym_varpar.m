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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
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
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:24
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
	% StartTime: 2022-01-23 09:33:25
	% EndTime: 2022-01-23 09:33:25
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
	t139 = t150 * qJD(1) + t147 * qJD(3);
	t138 = t148 * qJD(1) + t149 * qJD(3);
	t137 = t149 * qJD(1) + t148 * qJD(3);
	t136 = t147 * qJD(1) + t150 * qJD(3);
	t1 = [t139 * r_i_i_C(1) + t138 * r_i_i_C(2) + t146 * qJD(2) + (-qJ(2) * t144 + t151 * t146) * qJD(1), qJD(1) * t146, t136 * r_i_i_C(1) + t137 * r_i_i_C(2), 0, 0; -t137 * r_i_i_C(1) + t136 * r_i_i_C(2) + t144 * qJD(2) + (qJ(2) * t146 + t151 * t144) * qJD(1), qJD(1) * t144, -t138 * r_i_i_C(1) + t139 * r_i_i_C(2), 0, 0; 0, 0, (-r_i_i_C(1) * t145 + r_i_i_C(2) * t143) * t141 * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:24
	% EndTime: 2022-01-23 09:33:25
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (155->40), mult. (215->63), div. (0->0), fcn. (180->8), ass. (0->35)
	t198 = pkin(3) * qJD(3);
	t177 = cos(pkin(8));
	t179 = sin(qJ(1));
	t197 = t177 * t179;
	t181 = cos(qJ(1));
	t196 = t177 * t181;
	t178 = sin(qJ(3));
	t195 = t178 * t179;
	t194 = t178 * t181;
	t180 = cos(qJ(3));
	t193 = t179 * t180;
	t192 = t180 * t181;
	t174 = qJD(3) + qJD(4);
	t175 = qJ(3) + qJ(4);
	t172 = sin(t175);
	t173 = cos(t175);
	t182 = t172 * t197 + t173 * t181;
	t185 = -t172 * t179 - t173 * t196;
	t164 = t182 * qJD(1) + t185 * t174;
	t183 = t172 * t196 - t173 * t179;
	t184 = -t172 * t181 + t173 * t197;
	t165 = t184 * qJD(1) + t183 * t174;
	t191 = t164 * r_i_i_C(1) + t165 * r_i_i_C(2);
	t166 = t183 * qJD(1) + t184 * t174;
	t167 = t185 * qJD(1) + t182 * t174;
	t190 = -t166 * r_i_i_C(1) + t167 * r_i_i_C(2);
	t189 = r_i_i_C(1) * t173 * t174;
	t188 = t180 * t198;
	t176 = sin(pkin(8));
	t187 = -pkin(1) - (t180 * pkin(3) + pkin(2)) * t177 + (-r_i_i_C(3) - pkin(7) - pkin(6)) * t176;
	t186 = t177 * t178 * t198;
	t171 = t178 * pkin(3) + qJ(2);
	t170 = qJD(2) + t188;
	t169 = t176 * t174 * t172 * r_i_i_C(2);
	t1 = [t179 * t186 + t167 * r_i_i_C(1) + t166 * r_i_i_C(2) + t170 * t181 + (-t171 * t179 + t187 * t181) * qJD(1), qJD(1) * t181, ((-t177 * t192 - t195) * qJD(3) + (t177 * t195 + t192) * qJD(1)) * pkin(3) + t191, t191, 0; -t181 * t186 - t165 * r_i_i_C(1) + t164 * r_i_i_C(2) + t170 * t179 + (t171 * t181 + t187 * t179) * qJD(1), qJD(1) * t179, ((-t177 * t193 + t194) * qJD(3) + (-t177 * t194 + t193) * qJD(1)) * pkin(3) + t190, t190, 0; 0, 0, t169 + (-t188 - t189) * t176, -t176 * t189 + t169, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:33:25
	% EndTime: 2022-01-23 09:33:25
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (226->50), mult. (292->72), div. (0->0), fcn. (233->8), ass. (0->38)
	t182 = qJD(3) + qJD(4);
	t183 = qJ(3) + qJ(4);
	t178 = sin(t183);
	t179 = cos(t183);
	t187 = sin(qJ(1));
	t185 = cos(pkin(8));
	t189 = cos(qJ(1));
	t202 = t185 * t189;
	t192 = t178 * t202 - t179 * t187;
	t203 = t185 * t187;
	t193 = -t178 * t189 + t179 * t203;
	t169 = t192 * qJD(1) + t193 * t182;
	t207 = pkin(4) * t178;
	t206 = pkin(3) * qJD(3);
	t205 = t179 * t182;
	t184 = sin(pkin(8));
	t204 = t182 * t184;
	t186 = sin(qJ(3));
	t175 = t186 * pkin(3) + t207;
	t201 = qJ(2) + t175;
	t191 = t178 * t203 + t179 * t189;
	t194 = -t178 * t187 - t179 * t202;
	t167 = t191 * qJD(1) + t194 * t182;
	t168 = t193 * qJD(1) + t192 * t182;
	t200 = t167 * r_i_i_C(1) + t168 * r_i_i_C(2);
	t170 = t194 * qJD(1) + t191 * t182;
	t199 = -t169 * r_i_i_C(1) + t170 * r_i_i_C(2);
	t188 = cos(qJ(3));
	t176 = t188 * pkin(3) + pkin(4) * t179;
	t198 = qJD(1) * t187;
	t197 = qJD(1) * t189;
	t172 = pkin(4) * t205 + t188 * t206;
	t196 = qJD(2) + t172;
	t171 = -t182 * t207 - t186 * t206;
	t195 = qJD(5) * t184 + t171 * t185;
	t190 = -(pkin(2) + t176) * t185 - pkin(1) + (-r_i_i_C(3) - qJ(5) - pkin(7) - pkin(6)) * t184;
	t173 = t178 * r_i_i_C(2) * t204;
	t1 = [t170 * r_i_i_C(1) + t169 * r_i_i_C(2) + t196 * t189 - t195 * t187 + (-t201 * t187 + t190 * t189) * qJD(1), t197, -t172 * t202 + t187 * t171 + (t175 * t203 + t176 * t189) * qJD(1) + t200, t167 * pkin(4) + t200, -t184 * t198; -t168 * r_i_i_C(1) + t167 * r_i_i_C(2) + t195 * t189 + t196 * t187 + (t190 * t187 + t201 * t189) * qJD(1), t198, -t172 * t203 - t189 * t171 + (-t175 * t202 + t176 * t187) * qJD(1) + t199, -t169 * pkin(4) + t199, t184 * t197; 0, 0, t173 + (-r_i_i_C(1) * t205 - t172) * t184, t173 + (-pkin(4) - r_i_i_C(1)) * t179 * t204, 0;];
	JaD_transl = t1;
end