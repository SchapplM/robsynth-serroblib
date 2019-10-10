% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRPRR3
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:42
	% EndTime: 2019-10-10 00:49:42
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:42
	% EndTime: 2019-10-10 00:49:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:42
	% EndTime: 2019-10-10 00:49:42
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(10);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:42
	% EndTime: 2019-10-10 00:49:43
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (43->18), mult. (68->31), div. (0->0), fcn. (42->6), ass. (0->13)
	t24 = sin(qJ(3));
	t25 = cos(qJ(3));
	t26 = (r_i_i_C(1) * t24 + r_i_i_C(2) * t25) * qJD(3);
	t33 = pkin(7) + r_i_i_C(3);
	t32 = qJD(1) * t24;
	t31 = qJD(1) * t25;
	t30 = qJD(3) * t24;
	t29 = qJD(3) * t25;
	t27 = -r_i_i_C(1) * t25 + r_i_i_C(2) * t24 - pkin(2);
	t23 = qJ(1) + pkin(10);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21 * t26 + (-cos(qJ(1)) * pkin(1) - t33 * t21 + t27 * t22) * qJD(1), 0, (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0, 0; -t22 * t26 + (-sin(qJ(1)) * pkin(1) + t33 * t22 + t27 * t21) * qJD(1), 0, (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t26, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:44
	% EndTime: 2019-10-10 00:49:44
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (113->22), mult. (182->34), div. (0->0), fcn. (135->8), ass. (0->19)
	t167 = sin(qJ(3));
	t168 = cos(qJ(3));
	t165 = sin(pkin(11));
	t166 = cos(pkin(11));
	t176 = r_i_i_C(1) * t166 - r_i_i_C(2) * t165 + pkin(3);
	t180 = r_i_i_C(3) + qJ(4);
	t173 = t176 * t167 - t180 * t168;
	t182 = qJD(1) * t173;
	t181 = t173 * qJD(3) - t167 * qJD(4);
	t179 = qJD(1) * t167;
	t178 = qJD(3) * t168;
	t175 = t165 * r_i_i_C(1) + t166 * r_i_i_C(2) + pkin(7);
	t174 = -t180 * t167 - t176 * t168;
	t171 = -pkin(2) + t174;
	t170 = t174 * qJD(3) + qJD(4) * t168;
	t164 = qJ(1) + pkin(10);
	t163 = cos(t164);
	t162 = sin(t164);
	t1 = [t181 * t162 + (-cos(qJ(1)) * pkin(1) - t175 * t162 + t171 * t163) * qJD(1), 0, t162 * t182 + t170 * t163, -t162 * t179 + t163 * t178, 0, 0; -t181 * t163 + (-sin(qJ(1)) * pkin(1) + t175 * t163 + t171 * t162) * qJD(1), 0, t170 * t162 - t163 * t182, t162 * t178 + t163 * t179, 0, 0; 0, 0, -t181, qJD(3) * t167, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:44
	% EndTime: 2019-10-10 00:49:44
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (262->46), mult. (307->74), div. (0->0), fcn. (244->10), ass. (0->36)
	t216 = sin(qJ(3));
	t207 = cos(pkin(11)) * pkin(4) + pkin(3);
	t212 = pkin(11) + qJ(5);
	t208 = sin(t212);
	t210 = cos(t212);
	t222 = r_i_i_C(1) * t210 - r_i_i_C(2) * t208 + t207;
	t217 = cos(qJ(3));
	t239 = r_i_i_C(3) + pkin(8) + qJ(4);
	t240 = t239 * t217;
	t219 = -t222 * t216 + t240;
	t244 = qJD(1) * t219;
	t231 = t216 * qJD(4);
	t243 = (-t207 * t216 + t240) * qJD(3) + t231;
	t213 = qJ(1) + pkin(10);
	t211 = cos(t213);
	t237 = t210 * t211;
	t236 = qJD(1) * t216;
	t235 = qJD(3) * t216;
	t234 = qJD(3) * t217;
	t233 = qJD(5) * t216;
	t232 = qJD(5) * t217;
	t230 = sin(pkin(11)) * pkin(4) + pkin(7);
	t229 = t239 * t216;
	t226 = -qJD(1) + t232;
	t225 = qJD(1) * t217 - qJD(5);
	t224 = r_i_i_C(1) * t208 + r_i_i_C(2) * t210;
	t223 = t226 * t208;
	t221 = -t207 * t217 - pkin(2) - t229;
	t209 = sin(t213);
	t220 = t225 * t209 + t211 * t235;
	t218 = qJD(4) * t217 + t224 * t233 + (-t222 * t217 - t229) * qJD(3);
	t206 = -t225 * t237 + (t210 * t235 + t223) * t209;
	t205 = t226 * t210 * t209 + (-t209 * t235 + t225 * t211) * t208;
	t204 = t220 * t210 + t211 * t223;
	t203 = t220 * t208 - t226 * t237;
	t1 = [t206 * r_i_i_C(1) + t205 * r_i_i_C(2) - t243 * t209 + (-cos(qJ(1)) * pkin(1) - t230 * t209 + t221 * t211) * qJD(1), 0, -t209 * t244 + t218 * t211, -t209 * t236 + t211 * t234, t203 * r_i_i_C(1) + t204 * r_i_i_C(2), 0; -t204 * r_i_i_C(1) + t203 * r_i_i_C(2) + t243 * t211 + (-sin(qJ(1)) * pkin(1) + t230 * t211 + t221 * t209) * qJD(1), 0, t218 * t209 + t211 * t244, t209 * t234 + t211 * t236, -t205 * r_i_i_C(1) + t206 * r_i_i_C(2), 0; 0, 0, t219 * qJD(3) - t224 * t232 + t231, t235, (t208 * t233 - t210 * t234) * r_i_i_C(2) + (-t208 * t234 - t210 * t233) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:49:44
	% EndTime: 2019-10-10 00:49:44
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (530->64), mult. (453->93), div. (0->0), fcn. (363->12), ass. (0->56)
	t244 = sin(qJ(3));
	t241 = pkin(11) + qJ(5);
	t237 = cos(t241);
	t231 = pkin(5) * t237 + cos(pkin(11)) * pkin(4) + pkin(3);
	t239 = qJ(6) + t241;
	t233 = sin(t239);
	t278 = r_i_i_C(2) * t233;
	t234 = cos(t239);
	t279 = r_i_i_C(1) * t234;
	t251 = t231 - t278 + t279;
	t245 = cos(qJ(3));
	t275 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
	t283 = t275 * t245;
	t247 = -t251 * t244 + t283;
	t288 = qJD(1) * t247;
	t235 = sin(t241);
	t274 = pkin(5) * qJD(5);
	t264 = t235 * t274;
	t265 = t244 * qJD(4);
	t287 = (-t231 * t244 + t283) * qJD(3) - t245 * t264 + t265;
	t243 = qJ(1) + pkin(10);
	t238 = cos(t243);
	t242 = qJD(5) + qJD(6);
	t257 = t242 * t245 - qJD(1);
	t285 = t238 * t257;
	t268 = qJD(1) * t245;
	t256 = -t242 + t268;
	t236 = sin(t243);
	t267 = qJD(3) * t244;
	t262 = t236 * t267;
	t282 = t256 * t238 - t262;
	t277 = r_i_i_C(2) * t234;
	t254 = r_i_i_C(1) * t233 + t277;
	t281 = t254 * t242 + t264;
	t280 = pkin(5) * t235;
	t276 = pkin(7) + t280 + sin(pkin(11)) * pkin(4);
	t272 = t242 * t244;
	t261 = t238 * t267;
	t248 = t256 * t236 + t261;
	t226 = t248 * t233 - t234 * t285;
	t227 = t233 * t285 + t248 * t234;
	t271 = t226 * r_i_i_C(1) + t227 * r_i_i_C(2);
	t253 = t257 * t236;
	t228 = t282 * t233 + t234 * t253;
	t229 = t233 * t253 - t282 * t234;
	t270 = -t228 * r_i_i_C(1) + t229 * r_i_i_C(2);
	t269 = qJD(1) * t244;
	t266 = qJD(3) * t245;
	t263 = t237 * t274;
	t260 = t275 * t244;
	t255 = -qJD(5) + t268;
	t252 = (-qJD(5) * t245 + qJD(1)) * t237;
	t250 = -t231 * t245 - pkin(2) - t260;
	t246 = qJD(4) * t245 + t281 * t244 + (-t251 * t245 - t260) * qJD(3);
	t230 = t272 * t278;
	t1 = [t238 * t263 + t229 * r_i_i_C(1) + t228 * r_i_i_C(2) - t287 * t236 + (-cos(qJ(1)) * pkin(1) - t276 * t236 + t250 * t238) * qJD(1), 0, -t236 * t288 + t246 * t238, -t236 * t269 + t238 * t266, (t238 * t252 + (t255 * t236 + t261) * t235) * pkin(5) + t271, t271; t236 * t263 - t227 * r_i_i_C(1) + t226 * r_i_i_C(2) + t287 * t238 + (-sin(qJ(1)) * pkin(1) + t276 * t238 + t250 * t236) * qJD(1), 0, t246 * t236 + t238 * t288, t236 * t266 + t238 * t269, (t236 * t252 + (-t255 * t238 + t262) * t235) * pkin(5) + t270, t270; 0, 0, t247 * qJD(3) - t281 * t245 + t265, t267, t230 + (-t242 * t279 - t263) * t244 + (-t254 - t280) * t266, -t266 * t277 + t230 + (-t233 * t266 - t234 * t272) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end