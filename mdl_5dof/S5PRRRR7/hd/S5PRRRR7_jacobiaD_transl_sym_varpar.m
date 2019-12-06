% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRR7
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:13:42
	% EndTime: 2019-12-05 17:13:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:13:42
	% EndTime: 2019-12-05 17:13:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:13:42
	% EndTime: 2019-12-05 17:13:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->2), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->4)
	t12 = sin(qJ(2));
	t13 = cos(qJ(2));
	t14 = qJD(2) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, cos(pkin(9)) * t14, 0, 0, 0; 0, sin(pkin(9)) * t14, 0, 0, 0; 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:13:42
	% EndTime: 2019-12-05 17:13:42
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (31->19), mult. (116->43), div. (0->0), fcn. (87->6), ass. (0->17)
	t151 = pkin(6) + r_i_i_C(3);
	t138 = sin(qJ(3));
	t141 = cos(qJ(2));
	t150 = t138 * t141;
	t140 = cos(qJ(3));
	t149 = t140 * t141;
	t148 = qJD(2) * t141;
	t139 = sin(qJ(2));
	t147 = qJD(3) * t139;
	t146 = r_i_i_C(1) * t138 + r_i_i_C(2) * t140;
	t145 = -r_i_i_C(1) * t140 + r_i_i_C(2) * t138 - pkin(2);
	t144 = t146 * t139;
	t143 = qJD(2) * t144;
	t142 = qJD(3) * t144 + (-t151 * t139 + t145 * t141) * qJD(2);
	t137 = cos(pkin(9));
	t136 = sin(pkin(9));
	t1 = [0, t142 * t137, t137 * t143 + ((-t136 * t138 - t137 * t149) * r_i_i_C(1) + (-t136 * t140 + t137 * t150) * r_i_i_C(2)) * qJD(3), 0, 0; 0, t142 * t136, t136 * t143 + ((-t136 * t149 + t137 * t138) * r_i_i_C(1) + (t136 * t150 + t137 * t140) * r_i_i_C(2)) * qJD(3), 0, 0; 0, -t146 * t141 * qJD(3) + (t145 * t139 + t151 * t141) * qJD(2), (t138 * t147 - t140 * t148) * r_i_i_C(2) + (-t138 * t148 - t140 * t147) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:13:42
	% EndTime: 2019-12-05 17:13:42
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (131->33), mult. (206->62), div. (0->0), fcn. (158->8), ass. (0->33)
	t172 = qJD(3) + qJD(4);
	t176 = sin(qJ(3));
	t173 = qJ(3) + qJ(4);
	t170 = sin(t173);
	t171 = cos(t173);
	t199 = r_i_i_C(2) * t171;
	t185 = r_i_i_C(1) * t170 + t199;
	t197 = pkin(3) * qJD(3);
	t202 = t185 * t172 + t176 * t197;
	t201 = r_i_i_C(1) * t171;
	t200 = r_i_i_C(2) * t170;
	t198 = r_i_i_C(3) + pkin(7) + pkin(6);
	t177 = sin(qJ(2));
	t196 = t172 * t177;
	t179 = cos(qJ(2));
	t195 = t172 * t179;
	t178 = cos(qJ(3));
	t194 = t178 * t179;
	t174 = sin(pkin(9));
	t175 = cos(pkin(9));
	t191 = qJD(2) * t177;
	t182 = t172 * t175 + t174 * t191;
	t188 = t174 * t195;
	t193 = (t182 * t170 - t171 * t188) * r_i_i_C(1) + (t170 * t188 + t182 * t171) * r_i_i_C(2);
	t183 = -t172 * t174 + t175 * t191;
	t187 = t175 * t195;
	t192 = (t183 * t170 - t171 * t187) * r_i_i_C(1) + (t170 * t187 + t183 * t171) * r_i_i_C(2);
	t190 = qJD(2) * t179;
	t186 = t176 * t191;
	t184 = -t178 * pkin(3) - pkin(2) + t200 - t201;
	t181 = t202 * t177 + (-t198 * t177 + t184 * t179) * qJD(2);
	t168 = t196 * t200;
	t1 = [0, t181 * t175, (t175 * t186 + (-t174 * t176 - t175 * t194) * qJD(3)) * pkin(3) + t192, t192, 0; 0, t181 * t174, (t174 * t186 + (-t174 * t194 + t175 * t176) * qJD(3)) * pkin(3) + t193, t193, 0; 0, -t202 * t179 + (t184 * t177 + t198 * t179) * qJD(2), t168 + (-t172 * t201 - t178 * t197) * t177 + (-pkin(3) * t176 - t185) * t190, -t190 * t199 + t168 + (-t170 * t190 - t171 * t196) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 17:13:42
	% EndTime: 2019-12-05 17:13:42
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (326->48), mult. (304->72), div. (0->0), fcn. (232->10), ass. (0->46)
	t185 = qJD(3) + qJD(4);
	t190 = sin(qJ(3));
	t214 = pkin(3) * qJD(3);
	t187 = qJ(3) + qJ(4);
	t182 = sin(t187);
	t219 = pkin(4) * t182;
	t174 = -t185 * t219 - t190 * t214;
	t181 = qJD(5) + t185;
	t184 = qJ(5) + t187;
	t179 = sin(t184);
	t180 = cos(t184);
	t216 = r_i_i_C(2) * t180;
	t199 = r_i_i_C(1) * t179 + t216;
	t220 = t199 * t181 - t174;
	t218 = r_i_i_C(1) * t180;
	t217 = r_i_i_C(2) * t179;
	t215 = r_i_i_C(3) + pkin(8) + pkin(7) + pkin(6);
	t191 = sin(qJ(2));
	t213 = t181 * t191;
	t193 = cos(qJ(2));
	t212 = t181 * t193;
	t183 = cos(t187);
	t211 = t183 * t185;
	t189 = cos(pkin(9));
	t188 = sin(pkin(9));
	t208 = qJD(2) * t191;
	t201 = t188 * t208;
	t195 = t181 * t189 + t201;
	t204 = t188 * t212;
	t210 = (t195 * t179 - t180 * t204) * r_i_i_C(1) + (t179 * t204 + t195 * t180) * r_i_i_C(2);
	t200 = t189 * t208;
	t196 = -t181 * t188 + t200;
	t203 = t189 * t212;
	t209 = (t196 * t179 - t180 * t203) * r_i_i_C(1) + (t179 * t203 + t196 * t180) * r_i_i_C(2);
	t207 = qJD(2) * t193;
	t206 = pkin(4) * t211;
	t205 = t181 * t218;
	t202 = t193 * t211;
	t192 = cos(qJ(3));
	t198 = -t192 * pkin(3) - pkin(4) * t183 - pkin(2) + t217 - t218;
	t175 = -t192 * t214 - t206;
	t178 = -t190 * pkin(3) - t219;
	t197 = t175 * t193 - t178 * t208;
	t194 = t220 * t191 + (-t215 * t191 + t198 * t193) * qJD(2);
	t176 = t213 * t217;
	t1 = [0, t194 * t189, t188 * t174 + t197 * t189 + t209, (-t189 * t202 + (-t185 * t188 + t200) * t182) * pkin(4) + t209, t209; 0, t194 * t188, -t189 * t174 + t197 * t188 + t210, (-t188 * t202 + (t185 * t189 + t201) * t182) * pkin(4) + t210, t210; 0, -t220 * t193 + (t198 * t191 + t215 * t193) * qJD(2), t176 + (t175 - t205) * t191 + (t178 - t199) * t207, t176 + (-t205 - t206) * t191 + (-t199 - t219) * t207, -t207 * t216 + t176 + (-t179 * t207 - t180 * t213) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end