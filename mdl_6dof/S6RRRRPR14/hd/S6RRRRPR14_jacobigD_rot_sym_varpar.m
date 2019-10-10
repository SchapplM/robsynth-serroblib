% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:54
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR14_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR14_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:13
	% EndTime: 2019-10-10 12:54:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (8->8), mult. (35->23), div. (0->0), fcn. (35->8), ass. (0->15)
	t139 = sin(pkin(6));
	t152 = t139 * cos(pkin(7));
	t142 = sin(qJ(2));
	t143 = sin(qJ(1));
	t151 = t142 * t143;
	t145 = cos(qJ(1));
	t150 = t142 * t145;
	t144 = cos(qJ(2));
	t149 = t143 * t144;
	t148 = t144 * t145;
	t147 = qJD(1) * t139;
	t138 = sin(pkin(7));
	t146 = qJD(2) * t138;
	t141 = cos(pkin(6));
	t1 = [0, t145 * t147, -(t141 * t151 - t148) * t146 + (-(-t141 * t148 + t151) * t138 + t145 * t152) * qJD(1), 0, 0, 0; 0, t143 * t147, -(-t141 * t150 - t149) * t146 + (-(-t141 * t149 - t150) * t138 + t143 * t152) * qJD(1), 0, 0, 0; 0, 0, t139 * t142 * t146, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:13
	% EndTime: 2019-10-10 12:54:13
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (38->24), mult. (136->54), div. (0->0), fcn. (140->10), ass. (0->31)
	t215 = sin(pkin(7));
	t216 = sin(pkin(6));
	t242 = t216 * t215;
	t217 = cos(pkin(7));
	t222 = cos(qJ(3));
	t241 = t217 * t222;
	t219 = sin(qJ(3));
	t223 = cos(qJ(2));
	t240 = t219 * t223;
	t220 = sin(qJ(2));
	t239 = t220 * t222;
	t221 = sin(qJ(1));
	t238 = t221 * t220;
	t237 = t221 * t223;
	t224 = cos(qJ(1));
	t236 = t224 * t220;
	t235 = t224 * t223;
	t234 = qJD(1) * t216;
	t233 = qJD(2) * t219;
	t232 = t221 * t242;
	t231 = t224 * t242;
	t230 = t221 * t234;
	t229 = t224 * t234;
	t218 = cos(pkin(6));
	t228 = t218 * t235 - t238;
	t227 = -t218 * t237 - t236;
	t226 = t218 * t236 + t237;
	t225 = t218 * t238 - t235;
	t214 = t227 * qJD(1) - t226 * qJD(2);
	t213 = -t228 * qJD(1) + t225 * qJD(2);
	t1 = [0, t229, -t213 * t215 + t217 * t229, -t213 * t241 + t227 * t233 + (-t226 * t219 - t222 * t231) * qJD(1) + (-t225 * t222 + (t227 * t217 + t232) * t219) * qJD(3), 0, 0; 0, t230, -t214 * t215 + t217 * t230, -t214 * t241 + t228 * t233 + (-t225 * t219 - t222 * t232) * qJD(1) + (t226 * t222 + (t228 * t217 - t231) * t219) * qJD(3), 0, 0; 0, 0, qJD(2) * t220 * t242, t218 * t215 * qJD(3) * t219 + ((t217 * t240 + t239) * qJD(3) + (t217 * t239 + t240) * qJD(2)) * t216, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:14
	% EndTime: 2019-10-10 12:54:14
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (38->24), mult. (136->54), div. (0->0), fcn. (140->10), ass. (0->31)
	t268 = sin(pkin(7));
	t269 = sin(pkin(6));
	t295 = t269 * t268;
	t270 = cos(pkin(7));
	t275 = cos(qJ(3));
	t294 = t270 * t275;
	t272 = sin(qJ(3));
	t276 = cos(qJ(2));
	t293 = t272 * t276;
	t273 = sin(qJ(2));
	t292 = t273 * t275;
	t274 = sin(qJ(1));
	t291 = t274 * t273;
	t290 = t274 * t276;
	t277 = cos(qJ(1));
	t289 = t277 * t273;
	t288 = t277 * t276;
	t287 = qJD(1) * t269;
	t286 = qJD(2) * t272;
	t285 = t274 * t295;
	t284 = t277 * t295;
	t283 = t274 * t287;
	t282 = t277 * t287;
	t271 = cos(pkin(6));
	t281 = t271 * t288 - t291;
	t280 = -t271 * t290 - t289;
	t279 = t271 * t289 + t290;
	t278 = t271 * t291 - t288;
	t267 = t280 * qJD(1) - t279 * qJD(2);
	t266 = -t281 * qJD(1) + t278 * qJD(2);
	t1 = [0, t282, -t266 * t268 + t270 * t282, -t266 * t294 + t280 * t286 + (-t279 * t272 - t275 * t284) * qJD(1) + (-t278 * t275 + (t280 * t270 + t285) * t272) * qJD(3), 0, 0; 0, t283, -t267 * t268 + t270 * t283, -t267 * t294 + t281 * t286 + (-t278 * t272 - t275 * t285) * qJD(1) + (t279 * t275 + (t281 * t270 - t284) * t272) * qJD(3), 0, 0; 0, 0, qJD(2) * t273 * t295, t271 * t268 * qJD(3) * t272 + ((t270 * t293 + t292) * qJD(3) + (t270 * t292 + t293) * qJD(2)) * t269, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:14
	% EndTime: 2019-10-10 12:54:15
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (100->49), mult. (332->104), div. (0->0), fcn. (355->12), ass. (0->46)
	t304 = sin(pkin(7));
	t312 = cos(qJ(4));
	t339 = t304 * t312;
	t305 = sin(pkin(6));
	t311 = sin(qJ(1));
	t338 = t305 * t311;
	t315 = cos(qJ(1));
	t337 = t305 * t315;
	t306 = cos(pkin(7));
	t309 = sin(qJ(3));
	t336 = t306 * t309;
	t310 = sin(qJ(2));
	t335 = t309 * t310;
	t314 = cos(qJ(2));
	t334 = t309 * t314;
	t313 = cos(qJ(3));
	t333 = t310 * t313;
	t332 = t310 * t315;
	t331 = t311 * t310;
	t330 = t311 * t314;
	t329 = t313 * t314;
	t328 = t314 * t315;
	t327 = qJD(1) * t305;
	t308 = sin(qJ(4));
	t326 = qJD(3) * t308;
	t325 = t311 * t327;
	t324 = t315 * t327;
	t323 = t304 * t325;
	t322 = t304 * t324;
	t307 = cos(pkin(6));
	t300 = t307 * t328 - t331;
	t321 = t300 * t306 - t304 * t337;
	t302 = -t307 * t330 - t332;
	t320 = t302 * t306 + t304 * t338;
	t319 = t306 * t334 + t333;
	t301 = t307 * t332 + t330;
	t318 = t307 * t331 - t328;
	t317 = t301 * t313 + t321 * t309;
	t316 = t320 * t309 - t313 * t318;
	t299 = -t318 * qJD(1) + t300 * qJD(2);
	t298 = t302 * qJD(1) - t301 * qJD(2);
	t297 = -t301 * qJD(1) + t302 * qJD(2);
	t296 = -t300 * qJD(1) + t318 * qJD(2);
	t295 = -t298 * t304 + t306 * t325;
	t294 = -t296 * t304 + t306 * t324;
	t1 = [0, t324, t294, t297 * t309 + (-t296 * t306 - t322) * t313 + t316 * qJD(3), 0, (t296 * t336 + t297 * t313 + t309 * t322) * t308 - t294 * t312 + (t316 * t312 + (-t302 * t304 + t306 * t338) * t308) * qJD(4) + (t309 * t318 + t320 * t313) * t326; 0, t325, t295, t299 * t309 + (-t298 * t306 - t323) * t313 + t317 * qJD(3), 0, (t298 * t336 + t299 * t313 + t309 * t323) * t308 - t295 * t312 + (t317 * t312 + (-t300 * t304 - t306 * t337) * t308) * qJD(4) + (-t301 * t309 + t321 * t313) * t326; 0, 0, t305 * qJD(2) * t310 * t304, qJD(3) * t304 * t307 * t309 + (t319 * qJD(3) + (t306 * t333 + t334) * qJD(2)) * t305, 0, (t304 * t313 * t326 + (t306 * t308 + t309 * t339) * qJD(4)) * t307 + ((-t304 * t314 * t308 + t319 * t312) * qJD(4) + (t306 * t329 - t335) * t326 + ((-t306 * t335 + t329) * t308 - t310 * t339) * qJD(2)) * t305;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end