% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RPRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
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
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
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
	% StartTime: 2019-10-10 01:48:55
	% EndTime: 2019-10-10 01:48:55
	% DurationCPUTime: 0.26s
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
	t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t203 * t234 + (-cos(qJ(1)) * pkin(1) - pkin(7) * t203 + t213 * t204) * qJD(1), 0, -t203 * t235 + t210 * t204, r_i_i_C(1) * t199 + r_i_i_C(2) * t200, 0, 0; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t204 * t234 + (-sin(qJ(1)) * pkin(1) + pkin(7) * t204 + t213 * t203) * qJD(1), 0, t210 * t203 + t204 * t235, -r_i_i_C(1) * t201 + r_i_i_C(2) * t202, 0, 0; 0, 0, t211 * qJD(3) - t216 * t223, (t206 * t224 - t208 * t225) * r_i_i_C(2) + (-t206 * t225 - t208 * t224) * r_i_i_C(1), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:55
	% EndTime: 2019-10-10 01:48:55
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (381->55), mult. (424->84), div. (0->0), fcn. (334->10), ass. (0->51)
	t241 = sin(qJ(3));
	t242 = cos(qJ(4));
	t232 = t242 * pkin(4) + pkin(3);
	t239 = qJ(4) + qJ(5);
	t235 = sin(t239);
	t274 = r_i_i_C(2) * t235;
	t236 = cos(t239);
	t275 = r_i_i_C(1) * t236;
	t250 = t232 - t274 + t275;
	t243 = cos(qJ(3));
	t272 = r_i_i_C(3) + pkin(9) + pkin(8);
	t279 = t272 * t243;
	t246 = -t250 * t241 + t279;
	t285 = qJD(1) * t246;
	t240 = sin(qJ(4));
	t271 = pkin(4) * qJD(4);
	t263 = t240 * t271;
	t284 = (-t232 * t241 + t279) * qJD(3) - t243 * t263;
	t237 = qJD(4) + qJD(5);
	t266 = qJD(1) * t243;
	t255 = -t237 + t266;
	t282 = t236 * t255;
	t281 = t240 * (-qJD(4) + t266);
	t256 = t237 * t243 - qJD(1);
	t265 = qJD(3) * t241;
	t278 = -t235 * t265 + t256 * t236;
	t273 = r_i_i_C(2) * t236;
	t252 = r_i_i_C(1) * t235 + t273;
	t277 = t252 * t237 + t263;
	t276 = pkin(4) * t240;
	t269 = t237 * t241;
	t238 = qJ(1) + pkin(10);
	t233 = sin(t238);
	t234 = cos(t238);
	t251 = t255 * t235;
	t227 = t233 * t251 - t278 * t234;
	t248 = t256 * t235 + t236 * t265;
	t228 = t233 * t282 + t248 * t234;
	t268 = t227 * r_i_i_C(1) + t228 * r_i_i_C(2);
	t229 = t278 * t233 + t234 * t251;
	t230 = t248 * t233 - t234 * t282;
	t267 = -t229 * r_i_i_C(1) + t230 * r_i_i_C(2);
	t264 = qJD(3) * t243;
	t262 = t242 * t271;
	t261 = pkin(7) + t276;
	t259 = t272 * t241;
	t249 = -t232 * t243 - pkin(2) - t259;
	t247 = t240 * t265 + (-qJD(4) * t243 + qJD(1)) * t242;
	t245 = t277 * t241 + (-t250 * t243 - t259) * qJD(3);
	t231 = t269 * t274;
	t1 = [t234 * t262 + t230 * r_i_i_C(1) + t229 * r_i_i_C(2) - t284 * t233 + (-cos(qJ(1)) * pkin(1) - t261 * t233 + t249 * t234) * qJD(1), 0, -t233 * t285 + t245 * t234, (t233 * t281 + t247 * t234) * pkin(4) + t268, t268, 0; t233 * t262 - t228 * r_i_i_C(1) + t227 * r_i_i_C(2) + t284 * t234 + (-sin(qJ(1)) * pkin(1) + t261 * t234 + t249 * t233) * qJD(1), 0, t245 * t233 + t234 * t285, (t247 * t233 - t234 * t281) * pkin(4) + t267, t267, 0; 0, 0, t246 * qJD(3) - t277 * t243, t231 + (-t237 * t275 - t262) * t241 + (-t252 - t276) * t264, -t264 * t273 + t231 + (-t235 * t264 - t236 * t269) * r_i_i_C(1), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:56
	% EndTime: 2019-10-10 01:48:56
	% DurationCPUTime: 0.50s
	% Computational Cost: add. (739->72), mult. (742->106), div. (0->0), fcn. (622->10), ass. (0->59)
	t303 = qJ(4) + qJ(5);
	t300 = cos(t303);
	t347 = r_i_i_C(3) + qJ(6);
	t359 = t347 * t300;
	t301 = qJD(4) + qJD(5);
	t299 = sin(t303);
	t304 = sin(qJ(4));
	t346 = pkin(4) * qJD(4);
	t348 = r_i_i_C(2) + pkin(9) + pkin(8);
	t311 = qJD(3) * t348 + qJD(6) * t299 - t304 * t346;
	t350 = pkin(5) + r_i_i_C(1);
	t333 = t350 * t299;
	t358 = (t333 - t359) * t301 - t311;
	t306 = cos(qJ(4));
	t296 = pkin(4) * t306 + pkin(3);
	t307 = cos(qJ(3));
	t305 = sin(qJ(3));
	t337 = qJD(3) * t305;
	t357 = -t296 * t337 + t307 * t311;
	t317 = -t299 * t347 - t300 * t350;
	t314 = -t296 + t317;
	t355 = t314 * t305 + t348 * t307;
	t354 = (qJD(1) * t307 - qJD(4)) * t304;
	t349 = pkin(4) * t304;
	t302 = qJ(1) + pkin(10);
	t298 = cos(t302);
	t345 = t298 * t299;
	t344 = t299 * t301;
	t343 = t299 * t307;
	t342 = t300 * t301;
	t341 = t300 * t307;
	t340 = t301 * t305;
	t297 = sin(t302);
	t339 = qJD(1) * t297;
	t338 = qJD(1) * t298;
	t336 = qJD(3) * t307;
	t335 = t300 * qJD(6);
	t334 = t306 * t346;
	t332 = t298 * t341;
	t331 = t301 * t343;
	t330 = t300 * t340;
	t329 = t305 * t335 + t336 * t359;
	t328 = pkin(7) + t349;
	t327 = t300 * t339;
	t324 = t297 * t337;
	t321 = t334 - t335;
	t320 = t297 * t299 + t332;
	t319 = -t296 * t307 - t305 * t348 - pkin(2);
	t318 = t297 * t342 + t299 * t338;
	t277 = t337 * t345 - t301 * t332 - t297 * t344 + (t297 * t343 + t298 * t300) * qJD(1);
	t278 = t307 * t327 + (t300 * t337 + t331) * t298 - t318;
	t316 = qJD(6) * t320 + t277 * t350 - t278 * t347;
	t279 = -t298 * t344 - t299 * t324 + t307 * t318 - t327;
	t280 = qJD(1) * t320 - t297 * t331 - t298 * t342 - t300 * t324;
	t315 = -(-t297 * t341 + t345) * qJD(6) + t347 * t280 - t350 * t279;
	t313 = t304 * t337 + (-qJD(4) * t307 + qJD(1)) * t306;
	t312 = qJD(3) * t314;
	t309 = t305 * t358 + t307 * t312;
	t1 = [t321 * t298 - t350 * t280 - t347 * t279 - t357 * t297 + (-cos(qJ(1)) * pkin(1) - t328 * t297 + t319 * t298) * qJD(1), 0, t309 * t298 - t339 * t355, (t297 * t354 + t298 * t313) * pkin(4) + t316, t316, -t277; t321 * t297 - t350 * t278 - t347 * t277 + t357 * t298 + (-sin(qJ(1)) * pkin(1) + t328 * t298 + t319 * t297) * qJD(1), 0, t309 * t297 + t338 * t355, (t297 * t313 - t298 * t354) * pkin(4) + t315, t315, t279; 0, 0, t305 * t312 - t307 * t358, (-t333 - t349) * t336 + (t301 * t317 - t334) * t305 + t329, -t350 * t330 + (-t336 * t350 - t340 * t347) * t299 + t329, t299 * t336 + t330;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end