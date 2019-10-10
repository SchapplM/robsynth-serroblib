% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:26:58
	% EndTime: 2019-10-10 01:26:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:26:58
	% EndTime: 2019-10-10 01:26:59
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
	% StartTime: 2019-10-10 01:26:58
	% EndTime: 2019-10-10 01:26:58
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
	% StartTime: 2019-10-10 01:26:59
	% EndTime: 2019-10-10 01:26:59
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
	t23 = qJ(1) + pkin(10);
	t22 = cos(t23);
	t21 = sin(t23);
	t1 = [t21 * t26 + (-cos(qJ(1)) * pkin(1) - t33 * t21 + t27 * t22) * qJD(1), 0, (t21 * t31 + t22 * t30) * r_i_i_C(2) + (t21 * t32 - t22 * t29) * r_i_i_C(1), 0, 0, 0; -t22 * t26 + (-sin(qJ(1)) * pkin(1) + t33 * t22 + t27 * t21) * qJD(1), 0, (t21 * t30 - t22 * t31) * r_i_i_C(2) + (-t21 * t29 - t22 * t32) * r_i_i_C(1), 0, 0, 0; 0, 0, -t26, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:27:00
	% EndTime: 2019-10-10 01:27:00
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
	t205 = qJ(1) + pkin(10);
	t204 = cos(t205);
	t203 = sin(t205);
	t202 = t212 * t203 - t204 * t232;
	t201 = t229 * t203 + t204 * t214;
	t200 = t203 * t232 + t212 * t204;
	t199 = t203 * t214 - t229 * t204;
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t203 * t234 + (-cos(qJ(1)) * pkin(1) - pkin(7) * t203 + t213 * t204) * qJD(1), 0, -t203 * t235 + t210 * t204, t199 * r_i_i_C(1) + t200 * r_i_i_C(2), 0, 0; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t204 * t234 + (-sin(qJ(1)) * pkin(1) + pkin(7) * t204 + t213 * t203) * qJD(1), 0, t210 * t203 + t204 * t235, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2), 0, 0; 0, 0, t211 * qJD(3) - t216 * t223, (t206 * t224 - t208 * t225) * r_i_i_C(2) + (-t206 * t225 - t208 * t224) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:27:00
	% EndTime: 2019-10-10 01:27:01
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (318->54), mult. (514->84), div. (0->0), fcn. (429->8), ass. (0->41)
	t250 = sin(qJ(3));
	t252 = cos(qJ(3));
	t282 = pkin(8) + r_i_i_C(2);
	t284 = t282 * t252;
	t287 = (-pkin(3) * t250 + t284) * qJD(3);
	t249 = sin(qJ(4));
	t251 = cos(qJ(4));
	t279 = r_i_i_C(3) + qJ(5);
	t281 = r_i_i_C(1) + pkin(4);
	t283 = t281 * t249 - t279 * t251;
	t286 = t283 * qJD(4) - qJD(5) * t249;
	t258 = -t279 * t249 - t281 * t251;
	t255 = -pkin(3) + t258;
	t254 = t255 * t250 + t284;
	t248 = qJ(1) + pkin(10);
	t247 = cos(t248);
	t278 = t247 * t249;
	t277 = t249 * t252;
	t276 = t251 * t252;
	t246 = sin(t248);
	t275 = qJD(1) * t246;
	t274 = qJD(1) * t247;
	t273 = qJD(3) * t250;
	t272 = qJD(3) * t252;
	t271 = qJD(4) * t249;
	t270 = qJD(4) * t251;
	t268 = t282 * t250;
	t265 = t251 * t275;
	t264 = t246 * t273;
	t263 = t246 * t271;
	t262 = t247 * t270;
	t261 = t246 * t249 + t247 * t276;
	t260 = t246 * t277 + t247 * t251;
	t259 = -pkin(3) * t252 - pkin(2) - t268;
	t256 = t246 * t270 + t249 * t274;
	t253 = t286 * t250 + (t255 * t252 - t268) * qJD(3);
	t235 = t261 * qJD(1) - t251 * t264 - t252 * t263 - t262;
	t234 = -t247 * t271 - t249 * t264 + t256 * t252 - t265;
	t233 = t252 * t265 + (t251 * t273 + t252 * t271) * t247 - t256;
	t232 = t260 * qJD(1) - t252 * t262 + t273 * t278 - t263;
	t1 = [-t260 * qJD(5) - t281 * t235 - t279 * t234 - t246 * t287 + (-cos(qJ(1)) * pkin(1) - t246 * pkin(7) + t259 * t247) * qJD(1), 0, t253 * t247 - t254 * t275, t261 * qJD(5) + t281 * t232 - t279 * t233, -t232, 0; -(t246 * t251 - t247 * t277) * qJD(5) - t281 * t233 - t279 * t232 + t247 * t287 + (-sin(qJ(1)) * pkin(1) + t247 * pkin(7) + t259 * t246) * qJD(1), 0, t253 * t246 + t254 * t274, -(-t246 * t276 + t278) * qJD(5) + t279 * t235 - t281 * t234, t234, 0; 0, 0, t254 * qJD(3) - t286 * t252, -t283 * t272 + (t258 * qJD(4) + t251 * qJD(5)) * t250, t249 * t272 + t250 * t270, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:27:00
	% EndTime: 2019-10-10 01:27:01
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (712->91), mult. (1162->142), div. (0->0), fcn. (1084->10), ass. (0->58)
	t298 = sin(qJ(3));
	t301 = cos(qJ(3));
	t324 = -r_i_i_C(3) - pkin(9) + pkin(8);
	t343 = t324 * t301;
	t349 = (-t298 * pkin(3) + t343) * qJD(3);
	t297 = sin(qJ(4));
	t300 = cos(qJ(4));
	t299 = cos(qJ(6));
	t296 = sin(qJ(6));
	t336 = -pkin(4) - pkin(5);
	t316 = t296 * r_i_i_C(2) + t336;
	t307 = t299 * r_i_i_C(1) - t316;
	t319 = -t296 * r_i_i_C(1) - qJ(5);
	t308 = t299 * r_i_i_C(2) - t319;
	t337 = t297 * t308 + t300 * t307 + pkin(3);
	t348 = t337 * t298 - t343;
	t347 = (qJD(4) - qJD(6)) * t298;
	t295 = qJ(1) + pkin(10);
	t293 = sin(t295);
	t294 = cos(t295);
	t333 = t297 * t301;
	t277 = t293 * t333 + t294 * t300;
	t332 = t300 * t301;
	t334 = t294 * t297;
	t278 = t293 * t332 - t334;
	t313 = t277 * t296 + t278 * t299;
	t314 = t277 * t299 - t278 * t296;
	t345 = (r_i_i_C(1) * t313 + r_i_i_C(2) * t314) * qJD(6);
	t309 = t296 * t297 + t299 * t300;
	t310 = t296 * t300 - t297 * t299;
	t344 = -t297 * qJD(5) - t324 * qJD(3) + (r_i_i_C(1) * t310 + r_i_i_C(2) * t309) * qJD(6) + (t297 * t307 - t300 * t308) * qJD(4);
	t326 = qJD(4) * t300;
	t330 = qJD(1) * t294;
	t305 = t293 * t326 + t297 * t330;
	t329 = qJD(3) * t298;
	t322 = t293 * t329;
	t331 = qJD(1) * t293;
	t323 = t300 * t331;
	t327 = qJD(4) * t297;
	t275 = -t294 * t327 - t297 * t322 + t301 * t305 - t323;
	t272 = t275 * t299;
	t328 = qJD(3) * t301;
	t321 = t294 * t326;
	t320 = t301 * t327;
	t315 = -(-t309 * t347 + t310 * t328) * r_i_i_C(1) - (t309 * t328 + t310 * t347) * r_i_i_C(2);
	t279 = -t293 * t300 + t294 * t333;
	t280 = t293 * t297 + t294 * t332;
	t312 = t279 * t299 - t280 * t296;
	t311 = t279 * t296 + t280 * t299;
	t306 = -pkin(3) * t301 - t324 * t298 - pkin(2);
	t303 = qJD(3) * t337;
	t302 = t344 * t298 - t301 * t303;
	t276 = qJD(1) * t280 - t293 * t320 - t300 * t322 - t321;
	t274 = t301 * t323 + (t300 * t329 + t320) * t294 - t305;
	t273 = qJD(1) * t277 - t293 * t327 - t301 * t321 + t329 * t334;
	t269 = qJD(6) * t312 - t273 * t296 - t274 * t299;
	t268 = -qJD(6) * t311 - t273 * t299 + t274 * t296;
	t1 = [-t272 * r_i_i_C(2) - t277 * qJD(5) + t319 * t275 - t307 * t276 + (-r_i_i_C(1) * t314 + r_i_i_C(2) * t313) * qJD(6) - t293 * t349 + (-cos(qJ(1)) * pkin(1) - t293 * pkin(7) + t306 * t294) * qJD(1), 0, t302 * t294 + t348 * t331, t280 * qJD(5) - t308 * t274 + t307 * t273 + (r_i_i_C(1) * t311 + r_i_i_C(2) * t312) * qJD(6), -t273, t268 * r_i_i_C(1) - t269 * r_i_i_C(2); t269 * r_i_i_C(1) + t268 * r_i_i_C(2) - t273 * qJ(5) + t279 * qJD(5) + t336 * t274 + t294 * t349 + (-sin(qJ(1)) * pkin(1) + pkin(7) * t294 + t306 * t293) * qJD(1), 0, t302 * t293 - t348 * t330, -t272 * r_i_i_C(1) + t278 * qJD(5) + t316 * t275 + t308 * t276 + t345, t275, (-t276 * t296 + t272) * r_i_i_C(1) + (-t275 * t296 - t276 * t299) * r_i_i_C(2) - t345; 0, 0, -t298 * t303 - t344 * t301, (t300 * qJ(5) + t336 * t297) * t328 + (qJD(5) * t300 + (-t297 * qJ(5) + t336 * t300) * qJD(4)) * t298 - t315, t297 * t328 + t298 * t326, t315;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end