% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PRRPR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:31
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PRRPR6_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PRRPR6_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR6_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:48
	% EndTime: 2019-10-24 10:31:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:48
	% EndTime: 2019-10-24 10:31:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:48
	% EndTime: 2019-10-24 10:31:48
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
	% StartTime: 2019-10-24 10:31:49
	% EndTime: 2019-10-24 10:31:49
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
	% StartTime: 2019-10-24 10:31:50
	% EndTime: 2019-10-24 10:31:50
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (152->37), mult. (504->70), div. (0->0), fcn. (488->10), ass. (0->34)
	t256 = sin(pkin(9));
	t259 = cos(pkin(9));
	t262 = sin(qJ(2));
	t260 = cos(pkin(5));
	t264 = cos(qJ(2));
	t276 = t260 * t264;
	t285 = -t256 * t262 + t259 * t276;
	t277 = t260 * t262;
	t250 = t256 * t264 + t259 * t277;
	t263 = cos(qJ(3));
	t257 = sin(pkin(5));
	t261 = sin(qJ(3));
	t279 = t257 * t261;
	t284 = -t250 * t263 + t259 * t279;
	t255 = sin(pkin(10));
	t258 = cos(pkin(10));
	t272 = t258 * r_i_i_C(1) - t255 * r_i_i_C(2) + pkin(3);
	t282 = r_i_i_C(3) + qJ(4);
	t283 = t282 * t261 + t272 * t263 + pkin(2);
	t278 = t257 * t263;
	t273 = qJD(2) * t257 * t264;
	t271 = t255 * r_i_i_C(1) + t258 * r_i_i_C(2) + pkin(7);
	t268 = t256 * t277 - t259 * t264;
	t270 = t256 * t279 - t263 * t268;
	t269 = t256 * t276 + t259 * t262;
	t267 = t260 * t261 + t262 * t278;
	t266 = qJD(2) * t283;
	t265 = t261 * qJD(4) + (-t272 * t261 + t282 * t263) * qJD(3);
	t247 = t269 * qJD(2);
	t245 = t285 * qJD(2);
	t243 = t267 * qJD(3) + t261 * t273;
	t241 = t270 * qJD(3) - t247 * t261;
	t239 = -t284 * qJD(3) + t245 * t261;
	t1 = [0, -t271 * t247 - t265 * t269 + t268 * t266, t270 * qJD(4) + t282 * (-t247 * t263 + (t256 * t278 + t261 * t268) * qJD(3)) - t272 * t241, t241, 0; 0, t271 * t245 - t250 * t266 + t265 * t285, -t284 * qJD(4) + t282 * (t245 * t263 + (-t250 * t261 - t259 * t278) * qJD(3)) - t272 * t239, t239, 0; 0, (t265 * t264 + (-t283 * t262 + t271 * t264) * qJD(2)) * t257, t267 * qJD(4) + t282 * (t263 * t273 + (t260 * t263 - t262 * t279) * qJD(3)) - t272 * t243, t243, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:31:50
	% EndTime: 2019-10-24 10:31:50
	% DurationCPUTime: 0.38s
	% Computational Cost: add. (336->69), mult. (819->125), div. (0->0), fcn. (824->12), ass. (0->51)
	t305 = sin(qJ(3));
	t307 = cos(qJ(3));
	t300 = pkin(10) + qJ(5);
	t298 = sin(t300);
	t299 = cos(t300);
	t320 = t298 * r_i_i_C(1) + t299 * r_i_i_C(2);
	t315 = qJD(5) * t320;
	t321 = r_i_i_C(1) * t299 - r_i_i_C(2) * t298;
	t318 = cos(pkin(10)) * pkin(4) + pkin(3) + t321;
	t337 = r_i_i_C(3) + pkin(8) + qJ(4);
	t309 = (t318 * t305 - t337 * t307) * qJD(3) - t305 * qJD(4) + t307 * t315;
	t302 = sin(pkin(9));
	t306 = sin(qJ(2));
	t308 = cos(qJ(2));
	t335 = cos(pkin(9));
	t336 = cos(pkin(5));
	t319 = t336 * t335;
	t288 = t302 * t308 + t306 * t319;
	t303 = sin(pkin(5));
	t325 = t303 * t335;
	t278 = t288 * t307 - t305 * t325;
	t326 = t302 * t336;
	t290 = -t306 * t326 + t335 * t308;
	t333 = t303 * t305;
	t332 = t303 * t307;
	t331 = t303 * t308;
	t330 = qJD(2) * t306;
	t328 = t303 * t330;
	t327 = qJD(2) * t331;
	t317 = t308 * t319;
	t316 = -t290 * t305 + t302 * t332;
	t280 = t290 * t307 + t302 * t333;
	t314 = t321 * qJD(5);
	t313 = -t288 * t305 - t307 * t325;
	t312 = -t306 * t333 + t336 * t307;
	t292 = t336 * t305 + t306 * t332;
	t311 = sin(pkin(10)) * pkin(4) + pkin(7) + t320;
	t289 = t335 * t306 + t308 * t326;
	t310 = -t337 * t305 - t318 * t307 - pkin(2);
	t287 = t302 * t306 - t317;
	t286 = t290 * qJD(2);
	t285 = t289 * qJD(2);
	t284 = t288 * qJD(2);
	t283 = -qJD(2) * t317 + t302 * t330;
	t282 = t312 * qJD(3) + t307 * t327;
	t281 = t292 * qJD(3) + t305 * t327;
	t276 = t316 * qJD(3) - t285 * t307;
	t275 = t280 * qJD(3) - t285 * t305;
	t274 = t313 * qJD(3) - t283 * t307;
	t273 = t278 * qJD(3) - t283 * t305;
	t1 = [0, -t311 * t285 + t310 * t286 + t309 * t289 + t290 * t314, qJD(4) * t280 - t318 * t275 + t337 * t276 - t316 * t315, t275, (-t276 * t298 + t286 * t299) * r_i_i_C(1) + (-t276 * t299 - t286 * t298) * r_i_i_C(2) + ((-t280 * t299 - t289 * t298) * r_i_i_C(1) + (t280 * t298 - t289 * t299) * r_i_i_C(2)) * qJD(5); 0, -t311 * t283 + t310 * t284 + t309 * t287 + t288 * t314, qJD(4) * t278 - t318 * t273 + t337 * t274 - t313 * t315, t273, (-t274 * t298 + t284 * t299) * r_i_i_C(1) + (-t274 * t299 - t284 * t298) * r_i_i_C(2) + ((-t278 * t299 - t287 * t298) * r_i_i_C(1) + (t278 * t298 - t287 * t299) * r_i_i_C(2)) * qJD(5); 0, ((t310 * qJD(2) + t314) * t306 + (t311 * qJD(2) - t309) * t308) * t303, qJD(4) * t292 - t318 * t281 + t337 * t282 - t312 * t315, t281, (-t282 * t298 + t299 * t328) * r_i_i_C(1) + (-t282 * t299 - t298 * t328) * r_i_i_C(2) + ((-t292 * t299 + t298 * t331) * r_i_i_C(1) + (t292 * t298 + t299 * t331) * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end