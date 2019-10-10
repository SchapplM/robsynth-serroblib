% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (24->13), mult. (82->26), div. (0->0), fcn. (72->6), ass. (0->16)
	t135 = r_i_i_C(2) - pkin(2);
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
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:10
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-09 21:39:10
	% EndTime: 2019-10-09 21:39:11
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (164->42), mult. (541->75), div. (0->0), fcn. (524->10), ass. (0->37)
	t263 = cos(pkin(6));
	t264 = sin(qJ(4));
	t266 = cos(qJ(4));
	t260 = sin(pkin(6));
	t267 = cos(qJ(2));
	t284 = t260 * t267;
	t289 = -t263 * t266 + t264 * t284;
	t258 = sin(pkin(11));
	t261 = cos(pkin(11));
	t274 = t261 * r_i_i_C(1) - t258 * r_i_i_C(2) + pkin(4);
	t287 = r_i_i_C(3) + qJ(5);
	t288 = t274 * t264 - t287 * t266 + qJ(3);
	t286 = t260 * t264;
	t285 = t260 * t266;
	t265 = sin(qJ(2));
	t283 = t263 * t265;
	t281 = t263 * t267;
	t280 = qJD(2) * t265;
	t279 = qJD(2) * t267;
	t259 = sin(pkin(10));
	t277 = t259 * t280;
	t276 = t260 * t280;
	t262 = cos(pkin(10));
	t275 = t262 * t279;
	t273 = -t258 * r_i_i_C(1) - t261 * r_i_i_C(2) - pkin(2) - pkin(8);
	t250 = t259 * t265 - t262 * t281;
	t272 = -t250 * t264 + t262 * t285;
	t252 = t259 * t281 + t262 * t265;
	t271 = t252 * t264 + t259 * t285;
	t270 = t259 * t267 + t262 * t283;
	t268 = -t266 * qJD(5) + qJD(3) + (t287 * t264 + t274 * t266) * qJD(4);
	t249 = -t263 * t277 + t275;
	t247 = t270 * qJD(2);
	t244 = -t289 * qJD(4) - t266 * t276;
	t242 = t272 * qJD(4) + t247 * t266;
	t240 = t271 * qJD(4) - t249 * t266;
	t1 = [0, t273 * t249 - t288 * t252 * qJD(2) + t268 * (-t259 * t283 + t262 * t267), t249, t271 * qJD(5) - t287 * (-t249 * t264 + (-t252 * t266 + t259 * t286) * qJD(4)) - t274 * t240, t240, 0; 0, t273 * t247 - t288 * (-t263 * t275 + t277) + t268 * t270, t247, -t272 * qJD(5) + t287 * (t247 * t264 + (t250 * t266 + t262 * t286) * qJD(4)) + t274 * t242, -t242, 0; 0, (t288 * t279 + (t273 * qJD(2) + t268) * t265) * t260, t276, -t289 * qJD(5) - t287 * (-t264 * t276 + (t263 * t264 + t266 * t284) * qJD(4)) - t274 * t244, t244, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:10
	% EndTime: 2019-10-09 21:39:11
	% DurationCPUTime: 0.55s
	% Computational Cost: add. (348->75), mult. (856->131), div. (0->0), fcn. (860->12), ass. (0->54)
	t306 = pkin(11) + qJ(6);
	t304 = sin(t306);
	t305 = cos(t306);
	t326 = r_i_i_C(1) * t305 - r_i_i_C(2) * t304;
	t321 = t326 * qJD(6);
	t313 = sin(qJ(4));
	t315 = cos(qJ(4));
	t324 = cos(pkin(11)) * pkin(5) + pkin(4) + t326;
	t341 = r_i_i_C(3) + pkin(9) + qJ(5);
	t342 = t324 * t313 - t341 * t315 + qJ(3);
	t309 = sin(pkin(6));
	t340 = t309 * t313;
	t314 = sin(qJ(2));
	t339 = t309 * t314;
	t338 = t309 * t315;
	t316 = cos(qJ(2));
	t337 = t309 * t316;
	t311 = cos(pkin(6));
	t336 = t311 * t314;
	t335 = t311 * t316;
	t334 = qJD(2) * t314;
	t333 = qJD(2) * t316;
	t332 = qJD(4) * t315;
	t308 = sin(pkin(10));
	t331 = t308 * t334;
	t330 = t309 * t333;
	t310 = cos(pkin(10));
	t329 = t310 * t333;
	t328 = qJD(4) * t340;
	t327 = t309 * t334;
	t325 = -r_i_i_C(1) * t304 - r_i_i_C(2) * t305;
	t291 = t308 * t314 - t310 * t335;
	t280 = -t291 * t313 + t310 * t338;
	t323 = t291 * t315 + t310 * t340;
	t293 = t308 * t335 + t310 * t314;
	t278 = t293 * t313 + t308 * t338;
	t292 = t308 * t316 + t310 * t336;
	t322 = t311 * t313 + t315 * t337;
	t320 = qJD(6) * t325;
	t319 = -pkin(5) * sin(pkin(11)) - pkin(2) - pkin(8) + t325;
	t317 = -t315 * qJD(5) + qJD(3) + t313 * t320 + (t341 * t313 + t324 * t315) * qJD(4);
	t296 = t311 * t315 - t313 * t337;
	t294 = -t308 * t336 + t310 * t316;
	t290 = -t311 * t331 + t329;
	t289 = t293 * qJD(2);
	t288 = t292 * qJD(2);
	t287 = -t311 * t329 + t331;
	t282 = t311 * t332 - t315 * t327 - t316 * t328;
	t281 = t322 * qJD(4) - t313 * t327;
	t276 = t280 * qJD(4) + t288 * t315;
	t275 = t323 * qJD(4) + t288 * t313;
	t274 = t278 * qJD(4) - t290 * t315;
	t273 = -t290 * t313 - t293 * t332 + t308 * t328;
	t1 = [0, -t289 * t342 + t319 * t290 - t293 * t321 + t317 * t294, t290, qJD(5) * t278 - t341 * t273 + (t293 * t315 - t308 * t340) * t320 - t324 * t274, t274, (t273 * t304 - t289 * t305) * r_i_i_C(1) + (t273 * t305 + t289 * t304) * r_i_i_C(2) + ((-t278 * t305 - t294 * t304) * r_i_i_C(1) + (t278 * t304 - t294 * t305) * r_i_i_C(2)) * qJD(6); 0, -t287 * t342 + t319 * t288 - t291 * t321 + t317 * t292, t288, -qJD(5) * t280 + t341 * t275 + t324 * t276 + t323 * t320, -t276, (-t275 * t304 - t287 * t305) * r_i_i_C(1) + (-t275 * t305 + t287 * t304) * r_i_i_C(2) + ((t280 * t305 - t292 * t304) * r_i_i_C(1) + (-t280 * t304 - t292 * t305) * r_i_i_C(2)) * qJD(6); 0, ((t342 * qJD(2) + t321) * t316 + (t319 * qJD(2) + t317) * t314) * t309, t327, qJD(5) * t296 - t341 * t281 - t324 * t282 - t322 * t320, t282, (t281 * t304 + t305 * t330) * r_i_i_C(1) + (t281 * t305 - t304 * t330) * r_i_i_C(2) + ((-t296 * t305 - t304 * t339) * r_i_i_C(1) + (t296 * t304 - t305 * t339) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end