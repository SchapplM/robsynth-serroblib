% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:50
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRP5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:50:08
	% EndTime: 2019-10-09 21:50:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:50:08
	% EndTime: 2019-10-09 21:50:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:50:08
	% EndTime: 2019-10-09 21:50:08
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
	% StartTime: 2019-10-09 21:50:09
	% EndTime: 2019-10-09 21:50:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (24->13), mult. (82->26), div. (0->0), fcn. (72->6), ass. (0->16)
	t135 = -pkin(2) + r_i_i_C(2);
	t134 = r_i_i_C(3) + qJ(3);
	t126 = cos(pkin(6));
	t127 = sin(qJ(2));
	t133 = t126 * t127;
	t128 = cos(qJ(2));
	t132 = t126 * t128;
	t131 = qJD(2) * t134;
	t123 = sin(pkin(10));
	t125 = cos(pkin(10));
	t130 = t123 * t128 + t125 * t133;
	t129 = t123 * t133 - t125 * t128;
	t124 = sin(pkin(6));
	t122 = t129 * qJD(2);
	t120 = t130 * qJD(2);
	t1 = [0, -t129 * qJD(3) - t135 * t122 - (t123 * t132 + t125 * t127) * t131, -t122, 0, 0, 0; 0, t130 * qJD(3) + t135 * t120 - (t123 * t127 - t125 * t132) * t131, t120, 0, 0, 0; 0, (qJD(3) * t127 + (t135 * t127 + t134 * t128) * qJD(2)) * t124, t124 * qJD(2) * t127, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:50:09
	% EndTime: 2019-10-09 21:50:09
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (67->29), mult. (228->61), div. (0->0), fcn. (208->8), ass. (0->28)
	t178 = sin(pkin(6));
	t181 = sin(qJ(4));
	t199 = t178 * t181;
	t183 = cos(qJ(4));
	t198 = t178 * t183;
	t184 = cos(qJ(2));
	t197 = t178 * t184;
	t179 = cos(pkin(10));
	t196 = t179 * t184;
	t180 = cos(pkin(6));
	t182 = sin(qJ(2));
	t195 = t180 * t182;
	t194 = t180 * t184;
	t193 = qJD(2) * t182;
	t192 = -pkin(2) - pkin(8) - r_i_i_C(3);
	t177 = sin(pkin(10));
	t191 = t177 * t193;
	t190 = t178 * t193;
	t189 = qJD(2) * t196;
	t188 = t183 * r_i_i_C(1) - t181 * r_i_i_C(2);
	t187 = t181 * r_i_i_C(1) + t183 * r_i_i_C(2) + qJ(3);
	t186 = t177 * t184 + t179 * t195;
	t173 = t177 * t194 + t179 * t182;
	t185 = t188 * qJD(4) + qJD(3);
	t171 = t177 * t182 - t179 * t194;
	t170 = -t180 * t191 + t189;
	t168 = t186 * qJD(2);
	t1 = [0, t185 * (-t177 * t195 + t196) + t192 * t170 - t187 * t173 * qJD(2), t170, t188 * t170 + ((-t173 * t181 - t177 * t198) * r_i_i_C(1) + (-t173 * t183 + t177 * t199) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t185 * t186 + t192 * t168 - t187 * (-t180 * t189 + t191), t168, t188 * t168 + ((-t171 * t181 + t179 * t198) * r_i_i_C(1) + (-t171 * t183 - t179 * t199) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t185 * t182 + (t192 * t182 + t187 * t184) * qJD(2)) * t178, t190, t188 * t190 + ((-t180 * t183 + t181 * t197) * r_i_i_C(1) + (t180 * t181 + t183 * t197) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:50:10
	% EndTime: 2019-10-09 21:50:10
	% DurationCPUTime: 0.49s
	% Computational Cost: add. (236->64), mult. (762->121), div. (0->0), fcn. (760->10), ass. (0->49)
	t293 = sin(pkin(10));
	t295 = cos(pkin(10));
	t299 = sin(qJ(2));
	t296 = cos(pkin(6));
	t302 = cos(qJ(2));
	t322 = t296 * t302;
	t285 = t293 * t322 + t295 * t299;
	t301 = cos(qJ(4));
	t294 = sin(pkin(6));
	t298 = sin(qJ(4));
	t327 = t294 * t298;
	t331 = -t285 * t301 + t293 * t327;
	t297 = sin(qJ(5));
	t300 = cos(qJ(5));
	t314 = r_i_i_C(1) * t300 - r_i_i_C(2) * t297;
	t306 = t314 * qJD(5);
	t312 = pkin(4) + t314;
	t329 = pkin(9) + r_i_i_C(3);
	t330 = t312 * t298 - t329 * t301 + qJ(3);
	t326 = t294 * t299;
	t325 = t294 * t301;
	t324 = t294 * t302;
	t323 = t296 * t299;
	t321 = qJD(2) * t299;
	t320 = qJD(2) * t302;
	t318 = t293 * t321;
	t317 = t294 * t320;
	t316 = t295 * t320;
	t315 = t294 * t321;
	t313 = -r_i_i_C(1) * t297 - r_i_i_C(2) * t300;
	t311 = -pkin(2) - pkin(8) + t313;
	t283 = t293 * t299 - t295 * t322;
	t310 = -t283 * t298 + t295 * t325;
	t309 = t283 * t301 + t295 * t327;
	t272 = t285 * t298 + t293 * t325;
	t284 = t293 * t302 + t295 * t323;
	t308 = t296 * t298 + t301 * t324;
	t307 = -t296 * t301 + t298 * t324;
	t305 = qJD(5) * t313;
	t303 = qJD(3) + t298 * t305 + (t329 * t298 + t312 * t301) * qJD(4);
	t286 = -t293 * t323 + t295 * t302;
	t282 = -t296 * t318 + t316;
	t281 = t285 * qJD(2);
	t280 = t284 * qJD(2);
	t279 = -t296 * t316 + t318;
	t275 = t308 * qJD(4) - t298 * t315;
	t269 = t309 * qJD(4) + t280 * t298;
	t267 = t331 * qJD(4) - t282 * t298;
	t1 = [0, -t281 * t330 + t311 * t282 - t285 * t306 + t303 * t286, t282, -t329 * t267 - t331 * t305 + t312 * (-t272 * qJD(4) + t282 * t301), (t267 * t297 - t281 * t300) * r_i_i_C(1) + (t267 * t300 + t281 * t297) * r_i_i_C(2) + ((-t272 * t300 - t286 * t297) * r_i_i_C(1) + (t272 * t297 - t286 * t300) * r_i_i_C(2)) * qJD(5), 0; 0, -t279 * t330 + t311 * t280 - t283 * t306 + t303 * t284, t280, t329 * t269 + t309 * t305 + t312 * (t310 * qJD(4) + t280 * t301), (-t269 * t297 - t279 * t300) * r_i_i_C(1) + (-t269 * t300 + t279 * t297) * r_i_i_C(2) + ((-t284 * t297 + t300 * t310) * r_i_i_C(1) + (-t284 * t300 - t297 * t310) * r_i_i_C(2)) * qJD(5), 0; 0, ((t330 * qJD(2) + t306) * t302 + (t311 * qJD(2) + t303) * t299) * t294, t315, -t329 * t275 - t308 * t305 + t312 * (t307 * qJD(4) + t301 * t315), (t275 * t297 + t300 * t317) * r_i_i_C(1) + (t275 * t300 - t297 * t317) * r_i_i_C(2) + ((-t297 * t326 + t300 * t307) * r_i_i_C(1) + (-t297 * t307 - t300 * t326) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:50:10
	% EndTime: 2019-10-09 21:50:10
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (331->78), mult. (1015->134), div. (0->0), fcn. (1016->10), ass. (0->61)
	t311 = sin(qJ(5));
	t314 = cos(qJ(5));
	t350 = pkin(5) + r_i_i_C(1);
	t353 = -t311 * r_i_i_C(2) + t350 * t314;
	t321 = t353 * qJD(5);
	t312 = sin(qJ(4));
	t315 = cos(qJ(4));
	t327 = pkin(4) + t353;
	t348 = r_i_i_C(3) + qJ(6) + pkin(9);
	t351 = t327 * t312 - t348 * t315 + qJ(3);
	t347 = cos(pkin(6));
	t307 = sin(pkin(10));
	t309 = cos(pkin(10));
	t313 = sin(qJ(2));
	t316 = cos(qJ(2));
	t335 = t316 * t347;
	t296 = t307 * t335 + t309 * t313;
	t346 = t296 * t315;
	t308 = sin(pkin(6));
	t345 = t308 * t312;
	t344 = t308 * t313;
	t343 = t308 * t315;
	t342 = t308 * t316;
	t341 = qJD(2) * t313;
	t340 = qJD(2) * t316;
	t339 = t308 * t340;
	t338 = qJD(4) * t345;
	t337 = t308 * t341;
	t336 = t313 * t347;
	t334 = t347 * t315;
	t333 = t307 * t336;
	t332 = t309 * t335;
	t293 = -qJD(2) * t333 + t309 * t340;
	t276 = -qJD(4) * t346 - t293 * t312 + t307 * t338;
	t292 = t296 * qJD(2);
	t331 = t276 * t311 - t292 * t314;
	t295 = t307 * t316 + t309 * t336;
	t291 = t295 * qJD(2);
	t294 = t307 * t313 - t332;
	t325 = t294 * t315 + t309 * t345;
	t278 = t325 * qJD(4) + t291 * t312;
	t290 = -qJD(2) * t332 + t307 * t341;
	t330 = -t278 * t311 - t290 * t314;
	t281 = t296 * t312 + t307 * t343;
	t297 = t309 * t316 - t333;
	t329 = -t281 * t314 - t297 * t311;
	t283 = -t294 * t312 + t309 * t343;
	t328 = t283 * t314 - t295 * t311;
	t326 = -t314 * r_i_i_C(2) - t350 * t311;
	t299 = -t312 * t342 + t334;
	t324 = -t299 * t314 - t311 * t344;
	t323 = -pkin(2) - pkin(8) + t326;
	t319 = t347 * t312 + t315 * t342;
	t284 = t319 * qJD(4) - t312 * t337;
	t322 = t284 * t311 + t314 * t339;
	t320 = qJD(5) * t326;
	t317 = -t315 * qJD(6) + qJD(3) + t312 * t320 + (t348 * t312 + t327 * t315) * qJD(4);
	t285 = qJD(4) * t334 - t315 * t337 - t316 * t338;
	t279 = t283 * qJD(4) + t291 * t315;
	t277 = t281 * qJD(4) - t293 * t315;
	t1 = [0, -t292 * t351 + t323 * t293 - t296 * t321 + t317 * t297, t293, qJD(6) * t281 - t348 * t276 - t327 * t277 + (-t307 * t345 + t346) * t320, t331 * r_i_i_C(1) + (t276 * t314 + t292 * t311) * r_i_i_C(2) + (t329 * r_i_i_C(1) + (t281 * t311 - t297 * t314) * r_i_i_C(2)) * qJD(5) + (t329 * qJD(5) + t331) * pkin(5), t277; 0, -t290 * t351 + t323 * t291 - t294 * t321 + t317 * t295, t291, -qJD(6) * t283 + t348 * t278 + t327 * t279 + t325 * t320, t330 * r_i_i_C(1) + (-t278 * t314 + t290 * t311) * r_i_i_C(2) + (t328 * r_i_i_C(1) + (-t283 * t311 - t295 * t314) * r_i_i_C(2)) * qJD(5) + (t328 * qJD(5) + t330) * pkin(5), -t279; 0, ((t351 * qJD(2) + t321) * t316 + (t323 * qJD(2) + t317) * t313) * t308, t337, qJD(6) * t299 - t348 * t284 - t327 * t285 - t319 * t320, t322 * r_i_i_C(1) + (t284 * t314 - t311 * t339) * r_i_i_C(2) + (t324 * r_i_i_C(1) + (t299 * t311 - t314 * t344) * r_i_i_C(2)) * qJD(5) + (t324 * qJD(5) + t322) * pkin(5), t285;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end