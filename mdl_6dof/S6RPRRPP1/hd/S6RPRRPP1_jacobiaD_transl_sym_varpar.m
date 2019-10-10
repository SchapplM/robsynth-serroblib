% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:09
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
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
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->5), mult. (12->8), div. (0->0), fcn. (6->4), ass. (0->4)
	t32 = qJ(1) + pkin(9);
	t31 = cos(t32);
	t30 = sin(t32);
	t1 = [(-cos(qJ(1)) * pkin(1) - r_i_i_C(1) * t31 + r_i_i_C(2) * t30) * qJD(1), 0, 0, 0, 0, 0; (-sin(qJ(1)) * pkin(1) - r_i_i_C(1) * t30 - r_i_i_C(2) * t31) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:41
	% EndTime: 2019-10-10 01:09:41
	% DurationCPUTime: 0.10s
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
	t23 = qJ(1) + pkin(9);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21 * t26 + (-cos(qJ(1)) * pkin(1) - t33 * t21 + t27 * t22) * qJD(1), 0, (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0, 0; -t22 * t26 + (-sin(qJ(1)) * pkin(1) + t33 * t22 + t27 * t21) * qJD(1), 0, (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t26, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:43
	% EndTime: 2019-10-10 01:09:43
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (165->37), mult. (274->63), div. (0->0), fcn. (213->8), ass. (0->32)
	t207 = sin(qJ(3));
	t206 = sin(qJ(4));
	t208 = cos(qJ(4));
	t215 = r_i_i_C(1) * t208 - r_i_i_C(2) * t206 + pkin(3);
	t209 = cos(qJ(3));
	t228 = pkin(8) + r_i_i_C(3);
	t230 = t228 * t209;
	t211 = -t215 * t207 + t230;
	t235 = qJD(1) * t211;
	t234 = (-pkin(3) * t207 + t230) * qJD(3);
	t217 = qJD(1) * t209 - qJD(4);
	t232 = t208 * t217;
	t223 = qJD(4) * t209;
	t218 = -qJD(1) + t223;
	t226 = qJD(3) * t207;
	t229 = -t206 * t226 + t218 * t208;
	t225 = qJD(3) * t209;
	t224 = qJD(4) * t207;
	t222 = t228 * t207;
	t216 = r_i_i_C(1) * t206 + r_i_i_C(2) * t208;
	t214 = t217 * t206;
	t213 = -pkin(3) * t209 - pkin(2) - t222;
	t212 = t218 * t206 + t208 * t226;
	t210 = t216 * t224 + (-t215 * t209 - t222) * qJD(3);
	t205 = qJ(1) + pkin(9);
	t204 = cos(t205);
	t203 = sin(t205);
	t202 = t212 * t203 - t204 * t232;
	t201 = t229 * t203 + t204 * t214;
	t200 = t203 * t232 + t212 * t204;
	t199 = t203 * t214 - t229 * t204;
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t203 * t234 + (-cos(qJ(1)) * pkin(1) - pkin(7) * t203 + t213 * t204) * qJD(1), 0, -t203 * t235 + t210 * t204, r_i_i_C(1) * t199 + r_i_i_C(2) * t200, 0, 0; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t204 * t234 + (-sin(qJ(1)) * pkin(1) + pkin(7) * t204 + t213 * t203) * qJD(1), 0, t210 * t203 + t204 * t235, -r_i_i_C(1) * t201 + r_i_i_C(2) * t202, 0, 0; 0, 0, t211 * qJD(3) - t216 * t223, (t206 * t224 - t208 * t225) * r_i_i_C(2) + (-t206 * t225 - t208 * t224) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:43
	% EndTime: 2019-10-10 01:09:43
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (297->54), mult. (374->82), div. (0->0), fcn. (292->10), ass. (0->46)
	t218 = sin(qJ(3));
	t219 = cos(qJ(4));
	t249 = t219 * pkin(4);
	t209 = pkin(3) + t249;
	t214 = qJ(4) + pkin(10);
	t210 = sin(t214);
	t212 = cos(t214);
	t233 = r_i_i_C(1) * t212 - r_i_i_C(2) * t210;
	t229 = t209 + t233;
	t220 = cos(qJ(3));
	t248 = r_i_i_C(3) + qJ(5) + pkin(8);
	t251 = t248 * t220;
	t222 = -t229 * t218 + t251;
	t256 = qJD(1) * t222;
	t241 = t218 * qJD(5);
	t242 = qJD(4) * t220;
	t217 = sin(qJ(4));
	t250 = pkin(4) * t217;
	t255 = (-t209 * t218 + t251) * qJD(3) - t242 * t250 + t241;
	t215 = qJ(1) + pkin(9);
	t213 = cos(t215);
	t234 = qJD(1) * t220 - qJD(4);
	t231 = t234 * t213;
	t235 = -qJD(1) + t242;
	t252 = t235 * t212;
	t246 = qJD(1) * t218;
	t245 = qJD(3) * t218;
	t244 = qJD(3) * t220;
	t243 = qJD(4) * t218;
	t240 = qJD(4) * t249;
	t239 = pkin(7) + t250;
	t238 = t248 * t218;
	t232 = t235 * t210;
	t211 = sin(t215);
	t230 = t234 * t211;
	t228 = r_i_i_C(1) * t210 + r_i_i_C(2) * t212 + t250;
	t227 = -t209 * t220 - pkin(2) - t238;
	t225 = t228 * t220;
	t224 = t213 * t245 + t230;
	t223 = t217 * t245 - t235 * t219;
	t221 = qJD(5) * t220 + t228 * t243 + (-t229 * t220 - t238) * qJD(3);
	t208 = -t212 * t231 + (t212 * t245 + t232) * t211;
	t207 = t211 * t252 + (-t211 * t245 + t231) * t210;
	t206 = t224 * t212 + t213 * t232;
	t205 = t224 * t210 - t213 * t252;
	t1 = [t213 * t240 + t208 * r_i_i_C(1) + t207 * r_i_i_C(2) - t255 * t211 + (-cos(qJ(1)) * pkin(1) - t239 * t211 + t227 * t213) * qJD(1), 0, -t211 * t256 + t221 * t213, t205 * r_i_i_C(1) + t206 * r_i_i_C(2) + (t223 * t213 + t217 * t230) * pkin(4), -t211 * t246 + t213 * t244, 0; t211 * t240 - t206 * r_i_i_C(1) + t205 * r_i_i_C(2) + t255 * t213 + (-sin(qJ(1)) * pkin(1) + t239 * t213 + t227 * t211) * qJD(1), 0, t221 * t211 + t213 * t256, -t207 * r_i_i_C(1) + t208 * r_i_i_C(2) + (t223 * t211 - t217 * t231) * pkin(4), t211 * t244 + t213 * t246, 0; 0, 0, t222 * qJD(3) - qJD(4) * t225 + t241, (-t233 - t249) * t243 - qJD(3) * t225, t245, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:09:43
	% EndTime: 2019-10-10 01:09:44
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (534->71), mult. (614->104), div. (0->0), fcn. (508->10), ass. (0->55)
	t269 = cos(qJ(4));
	t308 = t269 * pkin(4);
	t259 = pkin(3) + t308;
	t267 = sin(qJ(4));
	t268 = sin(qJ(3));
	t270 = cos(qJ(3));
	t264 = qJ(4) + pkin(10);
	t260 = sin(t264);
	t307 = r_i_i_C(2) + qJ(5) + pkin(8);
	t279 = t307 * qJD(3) + t260 * qJD(6);
	t305 = pkin(4) * qJD(4);
	t319 = (-t267 * t305 + t279) * t270 - (qJD(3) * t259 - qJD(5)) * t268;
	t262 = cos(t264);
	t306 = r_i_i_C(3) + qJ(6);
	t309 = pkin(4) * t267;
	t310 = r_i_i_C(1) + pkin(5);
	t311 = t310 * t260 - t306 * t262 + t309;
	t317 = t311 * qJD(4) - t279;
	t277 = -t306 * t260 - t310 * t262;
	t275 = -t259 + t277;
	t316 = t275 * t268 + t307 * t270;
	t298 = qJD(1) * t270;
	t315 = t267 * (-qJD(4) + t298);
	t265 = qJ(1) + pkin(9);
	t261 = sin(t265);
	t304 = t261 * t260;
	t303 = t261 * t270;
	t263 = cos(t265);
	t302 = t263 * t262;
	t301 = qJD(1) * t261;
	t300 = qJD(1) * t263;
	t299 = qJD(1) * t268;
	t297 = qJD(3) * t268;
	t296 = qJD(3) * t270;
	t295 = qJD(4) * t262;
	t294 = qJD(4) * t263;
	t292 = t262 * qJD(6);
	t290 = pkin(7) + t309;
	t289 = t261 * t297;
	t288 = t263 * t297;
	t287 = qJD(4) * t304;
	t286 = t260 * t294;
	t285 = t262 * t294;
	t282 = t269 * t305 - t292;
	t281 = t270 * t302 + t304;
	t278 = -t259 * t270 - t307 * t268 - pkin(2);
	t276 = t260 * t300 + t261 * t295;
	t274 = t267 * t297 + (-qJD(4) * t270 + qJD(1)) * t269;
	t272 = t275 * qJD(3) + qJD(5);
	t271 = t317 * t268 + t272 * t270;
	t248 = t281 * qJD(1) - t262 * t289 - t270 * t287 - t285;
	t247 = -t260 * t289 - t262 * t301 + t276 * t270 - t286;
	t246 = t270 * t286 + (t261 * t298 + t288) * t262 - t276;
	t245 = t260 * t288 - t270 * t285 - t287 + (t260 * t303 + t302) * qJD(1);
	t1 = [t282 * t263 - t310 * t248 - t306 * t247 - t319 * t261 + (-cos(qJ(1)) * pkin(1) - t290 * t261 + t278 * t263) * qJD(1), 0, t271 * t263 - t316 * t301, t281 * qJD(6) - t306 * t246 + t310 * t245 + (t261 * t315 + t263 * t274) * pkin(4), -t261 * t299 + t263 * t296, -t245; t282 * t261 - t310 * t246 - t306 * t245 + t319 * t263 + (-sin(qJ(1)) * pkin(1) + t290 * t263 + t278 * t261) * qJD(1), 0, t271 * t261 + t316 * t300, -(t263 * t260 - t262 * t303) * qJD(6) + t306 * t248 - t310 * t247 + (t261 * t274 - t263 * t315) * pkin(4), t261 * t296 + t263 * t299, t247; 0, 0, t272 * t268 - t317 * t270, -t311 * t296 + (t292 + (t277 - t308) * qJD(4)) * t268, t297, t260 * t296 + t268 * t295;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end