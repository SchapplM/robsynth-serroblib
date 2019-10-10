% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:14
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:35
	% EndTime: 2019-10-09 22:14:35
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(10));
	t50 = sin(pkin(10));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:36
	% EndTime: 2019-10-09 22:14:36
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(10));
	t182 = cos(pkin(10));
	t185 = sin(qJ(2));
	t183 = cos(pkin(6));
	t187 = cos(qJ(2));
	t195 = t183 * t187;
	t201 = -t180 * t185 + t182 * t195;
	t200 = pkin(8) + r_i_i_C(3);
	t181 = sin(pkin(6));
	t184 = sin(qJ(3));
	t198 = t181 * t184;
	t186 = cos(qJ(3));
	t197 = t181 * t186;
	t196 = t183 * t185;
	t193 = t184 * r_i_i_C(1) + t186 * r_i_i_C(2);
	t192 = t186 * r_i_i_C(1) - t184 * r_i_i_C(2) + pkin(2);
	t176 = t180 * t187 + t182 * t196;
	t191 = t180 * t195 + t182 * t185;
	t190 = t180 * t196 - t182 * t187;
	t189 = qJD(3) * t193;
	t188 = qJD(2) * t192;
	t173 = t191 * qJD(2);
	t171 = t201 * qJD(2);
	t1 = [0, -t200 * t173 + t190 * t188 + t191 * t189, t193 * t173 + ((-t180 * t198 + t186 * t190) * r_i_i_C(1) + (-t180 * t197 - t184 * t190) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, t200 * t171 - t176 * t188 - t201 * t189, -t193 * t171 + ((-t176 * t186 + t182 * t198) * r_i_i_C(1) + (t176 * t184 + t182 * t197) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, (-t187 * t189 + (-t192 * t185 + t200 * t187) * qJD(2)) * t181, -t193 * t187 * t181 * qJD(2) + ((-t183 * t184 - t185 * t197) * r_i_i_C(1) + (-t183 * t186 + t185 * t198) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:36
	% EndTime: 2019-10-09 22:14:36
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (124->35), mult. (404->66), div. (0->0), fcn. (384->8), ass. (0->32)
	t221 = sin(pkin(10));
	t223 = cos(pkin(10));
	t226 = sin(qJ(2));
	t224 = cos(pkin(6));
	t228 = cos(qJ(2));
	t238 = t224 * t228;
	t247 = -t221 * t226 + t223 * t238;
	t225 = sin(qJ(3));
	t227 = cos(qJ(3));
	t243 = r_i_i_C(3) + qJ(4);
	t245 = pkin(3) - r_i_i_C(2);
	t246 = t243 * t225 + t245 * t227 + pkin(2);
	t244 = pkin(8) + r_i_i_C(1);
	t222 = sin(pkin(6));
	t241 = t222 * t225;
	t240 = t222 * t227;
	t239 = t224 * t226;
	t236 = qJD(2) * t222 * t228;
	t217 = t221 * t228 + t223 * t239;
	t235 = -t217 * t227 + t223 * t241;
	t232 = t221 * t239 - t223 * t228;
	t234 = t221 * t241 - t227 * t232;
	t233 = t221 * t238 + t223 * t226;
	t231 = t224 * t225 + t226 * t240;
	t230 = qJD(2) * t246;
	t229 = qJD(4) * t225 + (-t245 * t225 + t243 * t227) * qJD(3);
	t214 = t233 * qJD(2);
	t212 = t247 * qJD(2);
	t210 = t231 * qJD(3) + t225 * t236;
	t208 = t234 * qJD(3) - t214 * t225;
	t206 = -t235 * qJD(3) + t212 * t225;
	t1 = [0, -t244 * t214 - t229 * t233 + t232 * t230, t234 * qJD(4) + t243 * (-t214 * t227 + (t221 * t240 + t225 * t232) * qJD(3)) - t245 * t208, t208, 0, 0; 0, t244 * t212 - t217 * t230 + t229 * t247, -t235 * qJD(4) + t243 * (t212 * t227 + (-t217 * t225 - t223 * t240) * qJD(3)) - t245 * t206, t206, 0, 0; 0, (t229 * t228 + (-t246 * t226 + t244 * t228) * qJD(2)) * t222, t231 * qJD(4) + t243 * (t227 * t236 + (t224 * t227 - t226 * t241) * qJD(3)) - t245 * t210, t210, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:36
	% EndTime: 2019-10-09 22:14:37
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (203->43), mult. (659->75), div. (0->0), fcn. (644->10), ass. (0->41)
	t261 = sin(pkin(10));
	t264 = cos(pkin(10));
	t269 = cos(qJ(2));
	t265 = cos(pkin(6));
	t267 = sin(qJ(2));
	t287 = t265 * t267;
	t255 = t261 * t269 + t264 * t287;
	t268 = cos(qJ(3));
	t262 = sin(pkin(6));
	t266 = sin(qJ(3));
	t289 = t262 * t266;
	t292 = -t255 * t268 + t264 * t289;
	t260 = sin(pkin(11));
	t263 = cos(pkin(11));
	t280 = t260 * r_i_i_C(1) + t263 * r_i_i_C(2) + qJ(4);
	t284 = -r_i_i_C(3) - qJ(5) - pkin(3);
	t291 = t280 * t266 - t284 * t268 + pkin(2);
	t288 = t262 * t268;
	t286 = t265 * t269;
	t285 = qJD(2) * t267;
	t282 = t264 * t286;
	t281 = qJD(2) * t262 * t269;
	t279 = t263 * r_i_i_C(1) - t260 * r_i_i_C(2) + pkin(4) + pkin(8);
	t278 = t255 * t266 + t264 * t288;
	t274 = t261 * t287 - t264 * t269;
	t277 = t261 * t288 + t266 * t274;
	t276 = t261 * t289 - t268 * t274;
	t275 = t261 * t286 + t264 * t267;
	t273 = t265 * t266 + t267 * t288;
	t272 = -t265 * t268 + t267 * t289;
	t271 = qJD(2) * t291;
	t270 = t266 * qJD(4) + t268 * qJD(5) + (t284 * t266 + t280 * t268) * qJD(3);
	t252 = t275 * qJD(2);
	t250 = -qJD(2) * t282 + t261 * t285;
	t249 = -t272 * qJD(3) + t268 * t281;
	t248 = t273 * qJD(3) + t266 * t281;
	t247 = t277 * qJD(3) - t252 * t268;
	t246 = t276 * qJD(3) - t252 * t266;
	t245 = -t278 * qJD(3) - t250 * t268;
	t244 = -t292 * qJD(3) - t250 * t266;
	t1 = [0, -t279 * t252 - t270 * t275 + t274 * t271, t276 * qJD(4) + t277 * qJD(5) + t284 * t246 + t280 * t247, t246, t247, 0; 0, -t279 * t250 - t255 * t271 + t270 * (-t261 * t267 + t282), -t292 * qJD(4) - t278 * qJD(5) + t284 * t244 + t280 * t245, t244, t245, 0; 0, (-t291 * t285 + (t279 * qJD(2) + t270) * t269) * t262, t273 * qJD(4) - t272 * qJD(5) + t284 * t248 + t280 * t249, t248, t249, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:14:37
	% EndTime: 2019-10-09 22:14:37
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (399->71), mult. (1032->125), div. (0->0), fcn. (1036->12), ass. (0->51)
	t312 = sin(qJ(3));
	t314 = cos(qJ(3));
	t307 = pkin(11) + qJ(6);
	t305 = sin(t307);
	t306 = cos(t307);
	t329 = t306 * r_i_i_C(1) - t305 * r_i_i_C(2);
	t318 = t329 * qJD(6) + qJD(4);
	t328 = -t305 * r_i_i_C(1) - t306 * r_i_i_C(2);
	t319 = pkin(5) * sin(pkin(11)) + qJ(4) - t328;
	t337 = pkin(3) + r_i_i_C(3) + pkin(9) + qJ(5);
	t316 = (t337 * t312 - t319 * t314) * qJD(3) - t314 * qJD(5) - t318 * t312;
	t309 = sin(pkin(10));
	t313 = sin(qJ(2));
	t315 = cos(qJ(2));
	t344 = cos(pkin(10));
	t345 = cos(pkin(6));
	t327 = t345 * t344;
	t295 = t309 * t315 + t313 * t327;
	t310 = sin(pkin(6));
	t333 = t310 * t344;
	t347 = t295 * t314 - t312 * t333;
	t334 = t309 * t345;
	t297 = -t313 * t334 + t344 * t315;
	t342 = t310 * t312;
	t341 = t310 * t314;
	t340 = t310 * t315;
	t339 = qJD(2) * t313;
	t336 = t310 * t339;
	t335 = qJD(2) * t340;
	t326 = t315 * t327;
	t325 = -t297 * t312 + t309 * t341;
	t324 = t297 * t314 + t309 * t342;
	t323 = t328 * qJD(6);
	t322 = pkin(8) + cos(pkin(11)) * pkin(5) + pkin(4) + t329;
	t298 = t313 * t342 - t345 * t314;
	t321 = t345 * t312 + t313 * t341;
	t320 = -t295 * t312 - t314 * t333;
	t296 = t344 * t313 + t315 * t334;
	t317 = -t319 * t312 - t337 * t314 - pkin(2);
	t294 = t309 * t313 - t326;
	t293 = t297 * qJD(2);
	t292 = t296 * qJD(2);
	t291 = t295 * qJD(2);
	t290 = -qJD(2) * t326 + t309 * t339;
	t289 = -t298 * qJD(3) + t314 * t335;
	t288 = t321 * qJD(3) + t312 * t335;
	t283 = t325 * qJD(3) - t292 * t314;
	t282 = t324 * qJD(3) - t292 * t312;
	t281 = t320 * qJD(3) - t290 * t314;
	t280 = t347 * qJD(3) - t290 * t312;
	t1 = [0, -t322 * t292 + t317 * t293 + t316 * t296 + t297 * t323, qJD(5) * t325 - t337 * t282 + t319 * t283 + t318 * t324, t282, t283, (t282 * t306 - t293 * t305) * r_i_i_C(1) + (-t282 * t305 - t293 * t306) * r_i_i_C(2) + ((-t296 * t306 + t305 * t325) * r_i_i_C(1) + (t296 * t305 + t306 * t325) * r_i_i_C(2)) * qJD(6); 0, -t322 * t290 + t317 * t291 + t316 * t294 + t295 * t323, qJD(5) * t320 - t337 * t280 + t319 * t281 + t318 * t347, t280, t281, (t280 * t306 - t291 * t305) * r_i_i_C(1) + (-t280 * t305 - t291 * t306) * r_i_i_C(2) + ((-t294 * t306 + t305 * t320) * r_i_i_C(1) + (t294 * t305 + t306 * t320) * r_i_i_C(2)) * qJD(6); 0, ((t317 * qJD(2) + t323) * t313 + (t322 * qJD(2) - t316) * t315) * t310, -qJD(5) * t298 - t337 * t288 + t319 * t289 + t318 * t321, t288, t289, (t288 * t306 - t305 * t336) * r_i_i_C(1) + (-t288 * t305 - t306 * t336) * r_i_i_C(2) + ((-t298 * t305 + t306 * t340) * r_i_i_C(1) + (-t298 * t306 - t305 * t340) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end