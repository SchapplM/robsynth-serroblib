% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5RRPRR6
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRPRR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRPRR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
JaD_transl=NaN(3,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:06
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:05
	% EndTime: 2022-01-20 11:18:05
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (70->11), mult. (58->14), div. (0->0), fcn. (36->6), ass. (0->13)
	t37 = qJ(1) + qJ(2);
	t34 = sin(t37);
	t36 = qJD(1) + qJD(2);
	t46 = t36 * t34;
	t48 = qJ(3) + r_i_i_C(3);
	t47 = qJD(3) + r_i_i_C(2) * t36 * sin(pkin(9));
	t35 = cos(t37);
	t45 = t36 * t35;
	t44 = pkin(1) * qJD(1);
	t42 = -r_i_i_C(1) * cos(pkin(9)) - pkin(2);
	t41 = t34 * t47 + t42 * t46 + t45 * t48;
	t40 = -t48 * t46 + (t36 * t42 + t47) * t35;
	t1 = [-cos(qJ(1)) * t44 + t40, t40, t45, 0, 0; -sin(qJ(1)) * t44 + t41, t41, t46, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:06
	% EndTime: 2022-01-20 11:18:06
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (189->26), mult. (200->44), div. (0->0), fcn. (164->8), ass. (0->25)
	t189 = pkin(1) * qJD(1);
	t173 = qJ(1) + qJ(2);
	t170 = sin(t173);
	t172 = qJD(1) + qJD(2);
	t188 = t172 * t170;
	t171 = cos(t173);
	t187 = t172 * t171;
	t175 = cos(pkin(9));
	t176 = sin(qJ(4));
	t186 = t175 * t176;
	t177 = cos(qJ(4));
	t185 = t175 * t177;
	t184 = -t170 * t176 - t171 * t185;
	t183 = -t170 * t177 + t171 * t186;
	t182 = t170 * t185 - t171 * t176;
	t181 = t170 * t186 + t171 * t177;
	t174 = sin(pkin(9));
	t180 = -pkin(3) * t175 - pkin(2) + (-pkin(7) - r_i_i_C(3)) * t174;
	t163 = qJD(4) * t184 + t172 * t181;
	t164 = qJD(4) * t183 + t172 * t182;
	t179 = -t164 * r_i_i_C(1) + t163 * r_i_i_C(2) + qJ(3) * t187 + t170 * qJD(3) + t180 * t188;
	t165 = qJD(4) * t182 + t172 * t183;
	t166 = qJD(4) * t181 + t172 * t184;
	t178 = t165 * r_i_i_C(2) + t166 * r_i_i_C(1) + t171 * qJD(3) + (-qJ(3) * t170 + t171 * t180) * t172;
	t1 = [-cos(qJ(1)) * t189 + t178, t178, t187, r_i_i_C(1) * t163 + r_i_i_C(2) * t164, 0; -sin(qJ(1)) * t189 + t179, t179, t188, -r_i_i_C(1) * t165 + r_i_i_C(2) * t166, 0; 0, 0, 0, (-r_i_i_C(1) * t177 + r_i_i_C(2) * t176) * t174 * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-20 11:18:06
	% EndTime: 2022-01-20 11:18:06
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (389->42), mult. (325->64), div. (0->0), fcn. (262->10), ass. (0->40)
	t225 = sin(qJ(4));
	t251 = pkin(4) * t225 + qJ(3);
	t220 = qJD(1) + qJD(2);
	t223 = sin(pkin(9));
	t226 = cos(qJ(4));
	t247 = pkin(4) * qJD(4);
	t237 = t226 * t247;
	t250 = qJD(3) + t220 * t223 * (-pkin(8) - pkin(7)) + t237;
	t248 = pkin(1) * qJD(1);
	t221 = qJ(4) + qJ(5);
	t215 = sin(t221);
	t222 = qJ(1) + qJ(2);
	t218 = cos(t222);
	t246 = t215 * t218;
	t217 = cos(t221);
	t245 = t217 * t218;
	t216 = sin(t222);
	t244 = t220 * t216;
	t243 = t220 * t218;
	t224 = cos(pkin(9));
	t242 = t224 * t225;
	t241 = t224 * t226;
	t219 = qJD(4) + qJD(5);
	t235 = t220 * t224 - t219;
	t232 = t235 * t216;
	t234 = t219 * t224 - t220;
	t200 = t215 * t232 - t234 * t245;
	t201 = t217 * t232 + t234 * t246;
	t240 = t200 * r_i_i_C(1) + t201 * r_i_i_C(2);
	t231 = t234 * t216;
	t202 = t217 * t231 + t235 * t246;
	t203 = t215 * t231 - t235 * t245;
	t239 = -t202 * r_i_i_C(1) + t203 * r_i_i_C(2);
	t238 = r_i_i_C(1) * t217 * t219;
	t233 = t242 * t247;
	t230 = -r_i_i_C(3) * t223 - (t226 * pkin(4) + pkin(3)) * t224 - pkin(2);
	t229 = -t201 * r_i_i_C(1) + t200 * r_i_i_C(2) + t250 * t216 - t218 * t233 + t230 * t244 + t251 * t243;
	t228 = t203 * r_i_i_C(1) + t202 * r_i_i_C(2) + (-t251 * t220 + t233) * t216 + (t230 * t220 + t250) * t218;
	t207 = t223 * t219 * t215 * r_i_i_C(2);
	t1 = [-cos(qJ(1)) * t248 + t228, t228, t243, ((t216 * t242 + t218 * t226) * t220 + (-t216 * t225 - t218 * t241) * qJD(4)) * pkin(4) + t240, t240; -sin(qJ(1)) * t248 + t229, t229, t244, ((t216 * t226 - t218 * t242) * t220 + (-t216 * t241 + t218 * t225) * qJD(4)) * pkin(4) + t239, t239; 0, 0, 0, t207 + (-t237 - t238) * t223, -t223 * t238 + t207;];
	JaD_transl = t1;
end