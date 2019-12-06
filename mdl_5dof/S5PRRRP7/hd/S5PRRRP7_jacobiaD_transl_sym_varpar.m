% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRRP7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRRP7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRRP7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRP7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:06
	% EndTime: 2019-12-05 16:57:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:06
	% EndTime: 2019-12-05 16:57:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:07
	% EndTime: 2019-12-05 16:57:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(5));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(9));
	t50 = sin(pkin(9));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(5)) * qJD(2), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:07
	% EndTime: 2019-12-05 16:57:07
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(9));
	t182 = cos(pkin(9));
	t185 = sin(qJ(2));
	t183 = cos(pkin(5));
	t187 = cos(qJ(2));
	t195 = t183 * t187;
	t201 = -t180 * t185 + t182 * t195;
	t200 = pkin(7) + r_i_i_C(3);
	t181 = sin(pkin(5));
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
	t1 = [0, -t200 * t173 + t190 * t188 + t191 * t189, t193 * t173 + ((-t180 * t198 + t186 * t190) * r_i_i_C(1) + (-t180 * t197 - t184 * t190) * r_i_i_C(2)) * qJD(3), 0, 0; 0, t200 * t171 - t176 * t188 - t201 * t189, -t193 * t171 + ((-t176 * t186 + t182 * t198) * r_i_i_C(1) + (t176 * t184 + t182 * t197) * r_i_i_C(2)) * qJD(3), 0, 0; 0, (-t187 * t189 + (-t192 * t185 + t200 * t187) * qJD(2)) * t181, -t193 * t187 * t181 * qJD(2) + ((-t183 * t184 - t185 * t197) * r_i_i_C(1) + (-t183 * t186 + t185 * t198) * r_i_i_C(2)) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:08
	% EndTime: 2019-12-05 16:57:08
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (224->70), mult. (725->134), div. (0->0), fcn. (724->10), ass. (0->50)
	t291 = sin(qJ(3));
	t294 = cos(qJ(3));
	t290 = sin(qJ(4));
	t293 = cos(qJ(4));
	t305 = r_i_i_C(1) * t293 - r_i_i_C(2) * t290;
	t303 = pkin(3) + t305;
	t324 = pkin(8) + r_i_i_C(3);
	t326 = (t303 * t291 - t324 * t294) * qJD(3);
	t304 = t290 * r_i_i_C(1) + t293 * r_i_i_C(2);
	t321 = cos(pkin(5));
	t288 = sin(pkin(5));
	t320 = t288 * t291;
	t319 = t288 * t294;
	t295 = cos(qJ(2));
	t318 = t288 * t295;
	t292 = sin(qJ(2));
	t317 = qJD(2) * t292;
	t316 = qJD(2) * t295;
	t315 = qJD(4) * t290;
	t314 = qJD(4) * t293;
	t313 = qJD(4) * t294;
	t312 = t288 * t317;
	t311 = t288 * t316;
	t310 = t292 * t321;
	t309 = t295 * t321;
	t287 = sin(pkin(9));
	t307 = t287 * t310;
	t289 = cos(pkin(9));
	t306 = t289 * t309;
	t279 = t287 * t295 + t289 * t310;
	t302 = -t279 * t291 - t289 * t319;
	t301 = -t279 * t294 + t289 * t320;
	t281 = t289 * t295 - t307;
	t300 = -t281 * t291 + t287 * t319;
	t271 = t281 * t294 + t287 * t320;
	t299 = qJD(4) * t304;
	t280 = t287 * t309 + t289 * t292;
	t298 = -t292 * t320 + t321 * t294;
	t283 = t321 * t291 + t292 * t319;
	t297 = -t324 * t291 - t303 * t294 - pkin(2);
	t296 = t304 * t313 + t326;
	t278 = t287 * t292 - t306;
	t277 = -qJD(2) * t307 + t289 * t316;
	t276 = t280 * qJD(2);
	t275 = t279 * qJD(2);
	t274 = -qJD(2) * t306 + t287 * t317;
	t273 = t298 * qJD(3) + t294 * t311;
	t267 = t300 * qJD(3) - t276 * t294;
	t265 = t302 * qJD(3) - t274 * t294;
	t1 = [0, (-t276 * t290 + t281 * t314) * r_i_i_C(1) + (-t276 * t293 - t281 * t315) * r_i_i_C(2) - t276 * pkin(7) + t297 * t277 + t296 * t280, t324 * t267 - t300 * t299 + t303 * (-t271 * qJD(3) + t276 * t291), (-t267 * t290 + t277 * t293) * r_i_i_C(1) + (-t267 * t293 - t277 * t290) * r_i_i_C(2) + ((-t271 * t293 - t280 * t290) * r_i_i_C(1) + (t271 * t290 - t280 * t293) * r_i_i_C(2)) * qJD(4), 0; 0, (-t274 * t290 + t279 * t314) * r_i_i_C(1) + (-t274 * t293 - t279 * t315) * r_i_i_C(2) - t274 * pkin(7) + t297 * t275 + t296 * t278, t324 * t265 - t302 * t299 + t303 * (t301 * qJD(3) + t274 * t291), (-t265 * t290 + t275 * t293) * r_i_i_C(1) + (-t265 * t293 - t275 * t290) * r_i_i_C(2) + ((-t278 * t290 + t293 * t301) * r_i_i_C(1) + (-t278 * t293 - t290 * t301) * r_i_i_C(2)) * qJD(4), 0; 0, ((t297 * qJD(2) + t305 * qJD(4)) * t292 + (qJD(2) * pkin(7) - t326 + t304 * (qJD(2) - t313)) * t295) * t288, t324 * t273 - t298 * t299 + t303 * (-t283 * qJD(3) - t291 * t311), (-t273 * t290 + t293 * t312) * r_i_i_C(1) + (-t273 * t293 - t290 * t312) * r_i_i_C(2) + ((-t283 * t293 + t290 * t318) * r_i_i_C(1) + (t283 * t290 + t293 * t318) * r_i_i_C(2)) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:57:08
	% EndTime: 2019-12-05 16:57:08
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (319->80), mult. (978->143), div. (0->0), fcn. (980->10), ass. (0->59)
	t299 = sin(qJ(4));
	t302 = cos(qJ(4));
	t344 = pkin(4) + r_i_i_C(1);
	t348 = -t299 * r_i_i_C(2) + t344 * t302;
	t346 = t302 * r_i_i_C(2) + t344 * t299;
	t300 = sin(qJ(3));
	t303 = cos(qJ(3));
	t316 = pkin(3) + t348;
	t341 = r_i_i_C(3) + qJ(5) + pkin(8);
	t347 = -(t316 * t300 - t341 * t303) * qJD(3) + t300 * qJD(5);
	t296 = sin(pkin(9));
	t301 = sin(qJ(2));
	t304 = cos(qJ(2));
	t339 = cos(pkin(9));
	t340 = cos(pkin(5));
	t321 = t340 * t339;
	t286 = t296 * t304 + t301 * t321;
	t297 = sin(pkin(5));
	t326 = t297 * t339;
	t276 = t286 * t303 - t300 * t326;
	t327 = t296 * t340;
	t288 = -t301 * t327 + t339 * t304;
	t337 = t297 * t300;
	t336 = t297 * t303;
	t335 = t297 * t304;
	t334 = qJD(2) * t301;
	t333 = qJD(4) * t299;
	t332 = qJD(4) * t302;
	t331 = qJD(4) * t303;
	t329 = t297 * t334;
	t328 = qJD(2) * t335;
	t315 = t304 * t321;
	t281 = -qJD(2) * t315 + t296 * t334;
	t308 = -t286 * t300 - t303 * t326;
	t272 = t308 * qJD(3) - t281 * t303;
	t282 = t286 * qJD(2);
	t320 = -t272 * t299 + t282 * t302;
	t287 = t339 * t301 + t304 * t327;
	t283 = t287 * qJD(2);
	t314 = -t288 * t300 + t296 * t336;
	t274 = t314 * qJD(3) - t283 * t303;
	t284 = t288 * qJD(2);
	t319 = -t274 * t299 + t284 * t302;
	t285 = t296 * t301 - t315;
	t318 = -t276 * t302 - t285 * t299;
	t278 = t288 * t303 + t296 * t337;
	t317 = -t278 * t302 - t287 * t299;
	t290 = t340 * t300 + t301 * t336;
	t313 = -t290 * t302 + t299 * t335;
	t307 = -t301 * t337 + t340 * t303;
	t280 = t307 * qJD(3) + t303 * t328;
	t310 = -t280 * t299 + t302 * t329;
	t309 = qJD(4) * t346;
	t306 = -t341 * t300 - t316 * t303 - pkin(2);
	t305 = t346 * t331 - t347;
	t279 = t290 * qJD(3) + t300 * t328;
	t273 = t278 * qJD(3) - t283 * t300;
	t271 = t276 * qJD(3) - t281 * t300;
	t1 = [0, (-t283 * t302 - t288 * t333) * r_i_i_C(2) - t283 * pkin(7) + t306 * t284 + t305 * t287 + t344 * (-t283 * t299 + t288 * t332), t278 * qJD(5) - t316 * t273 + t341 * t274 - t314 * t309, t319 * r_i_i_C(1) + (-t274 * t302 - t284 * t299) * r_i_i_C(2) + (t317 * r_i_i_C(1) + (t278 * t299 - t287 * t302) * r_i_i_C(2)) * qJD(4) + (t317 * qJD(4) + t319) * pkin(4), t273; 0, (-t281 * t302 - t286 * t333) * r_i_i_C(2) - t281 * pkin(7) + t306 * t282 + t305 * t285 + t344 * (-t281 * t299 + t286 * t332), t276 * qJD(5) - t316 * t271 + t341 * t272 - t308 * t309, t320 * r_i_i_C(1) + (-t272 * t302 - t282 * t299) * r_i_i_C(2) + (t318 * r_i_i_C(1) + (t276 * t299 - t285 * t302) * r_i_i_C(2)) * qJD(4) + (t318 * qJD(4) + t320) * pkin(4), t271; 0, ((t306 * qJD(2) + t348 * qJD(4)) * t301 + (qJD(2) * pkin(7) + t346 * (qJD(2) - t331) + t347) * t304) * t297, t290 * qJD(5) - t316 * t279 + t341 * t280 - t307 * t309, t310 * r_i_i_C(1) + (-t280 * t302 - t299 * t329) * r_i_i_C(2) + (t313 * r_i_i_C(1) + (t290 * t299 + t302 * t335) * r_i_i_C(2)) * qJD(4) + (t313 * qJD(4) + t310) * pkin(4), t279;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end