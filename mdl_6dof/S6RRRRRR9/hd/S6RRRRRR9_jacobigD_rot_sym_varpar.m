% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:32
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR9_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR9_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:11
	% EndTime: 2019-10-10 13:32:11
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-10 13:32:12
	% EndTime: 2019-10-10 13:32:12
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-10 13:32:13
	% EndTime: 2019-10-10 13:32:13
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (100->49), mult. (332->104), div. (0->0), fcn. (355->12), ass. (0->46)
	t301 = sin(pkin(7));
	t309 = cos(qJ(4));
	t336 = t301 * t309;
	t302 = sin(pkin(6));
	t308 = sin(qJ(1));
	t335 = t302 * t308;
	t312 = cos(qJ(1));
	t334 = t302 * t312;
	t303 = cos(pkin(7));
	t306 = sin(qJ(3));
	t333 = t303 * t306;
	t307 = sin(qJ(2));
	t332 = t306 * t307;
	t311 = cos(qJ(2));
	t331 = t306 * t311;
	t310 = cos(qJ(3));
	t330 = t307 * t310;
	t329 = t308 * t307;
	t328 = t308 * t311;
	t327 = t310 * t311;
	t326 = t312 * t307;
	t325 = t312 * t311;
	t324 = qJD(1) * t302;
	t305 = sin(qJ(4));
	t323 = qJD(3) * t305;
	t322 = t308 * t324;
	t321 = t312 * t324;
	t320 = t301 * t322;
	t319 = t301 * t321;
	t304 = cos(pkin(6));
	t297 = t304 * t325 - t329;
	t318 = t297 * t303 - t301 * t334;
	t299 = -t304 * t328 - t326;
	t317 = t299 * t303 + t301 * t335;
	t316 = t303 * t331 + t330;
	t298 = t304 * t326 + t328;
	t315 = t304 * t329 - t325;
	t314 = t298 * t310 + t318 * t306;
	t313 = t317 * t306 - t310 * t315;
	t296 = -t315 * qJD(1) + t297 * qJD(2);
	t295 = t299 * qJD(1) - t298 * qJD(2);
	t294 = -t298 * qJD(1) + t299 * qJD(2);
	t293 = -t297 * qJD(1) + t315 * qJD(2);
	t292 = -t295 * t301 + t303 * t322;
	t291 = -t293 * t301 + t303 * t321;
	t1 = [0, t321, t291, t294 * t306 + (-t293 * t303 - t319) * t310 + t313 * qJD(3), (t293 * t333 + t294 * t310 + t306 * t319) * t305 - t291 * t309 + (t313 * t309 + (-t299 * t301 + t303 * t335) * t305) * qJD(4) + (t306 * t315 + t317 * t310) * t323, 0; 0, t322, t292, t296 * t306 + (-t295 * t303 - t320) * t310 + t314 * qJD(3), (t295 * t333 + t296 * t310 + t306 * t320) * t305 - t292 * t309 + (t314 * t309 + (-t297 * t301 - t303 * t334) * t305) * qJD(4) + (-t298 * t306 + t318 * t310) * t323, 0; 0, 0, t302 * qJD(2) * t307 * t301, t304 * t301 * qJD(3) * t306 + (t316 * qJD(3) + (t303 * t330 + t331) * qJD(2)) * t302, (t301 * t310 * t323 + (t303 * t305 + t306 * t336) * qJD(4)) * t304 + ((-t301 * t311 * t305 + t316 * t309) * qJD(4) + (t303 * t327 - t332) * t323 + ((-t303 * t332 + t327) * t305 - t307 * t336) * qJD(2)) * t302, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:32:13
	% EndTime: 2019-10-10 13:32:13
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (162->49), mult. (528->104), div. (0->0), fcn. (570->12), ass. (0->49)
	t327 = sin(pkin(7));
	t335 = cos(qJ(4));
	t362 = t327 * t335;
	t328 = sin(pkin(6));
	t334 = sin(qJ(1));
	t361 = t328 * t334;
	t338 = cos(qJ(1));
	t360 = t328 * t338;
	t329 = cos(pkin(7));
	t332 = sin(qJ(3));
	t359 = t329 * t332;
	t333 = sin(qJ(2));
	t358 = t332 * t333;
	t337 = cos(qJ(2));
	t357 = t332 * t337;
	t336 = cos(qJ(3));
	t356 = t333 * t336;
	t355 = t333 * t338;
	t354 = t334 * t333;
	t353 = t334 * t337;
	t352 = t336 * t337;
	t351 = t337 * t338;
	t350 = qJD(1) * t328;
	t331 = sin(qJ(4));
	t349 = qJD(3) * t331;
	t348 = t334 * t350;
	t347 = t338 * t350;
	t346 = t327 * t348;
	t345 = t327 * t347;
	t330 = cos(pkin(6));
	t323 = t330 * t351 - t354;
	t344 = t323 * t329 - t327 * t360;
	t325 = -t330 * t353 - t355;
	t343 = t325 * t329 + t327 * t361;
	t342 = t329 * t357 + t356;
	t324 = t330 * t355 + t353;
	t341 = t330 * t354 - t351;
	t340 = t324 * t336 + t344 * t332;
	t339 = t343 * t332 - t336 * t341;
	t322 = -t341 * qJD(1) + t323 * qJD(2);
	t321 = t325 * qJD(1) - t324 * qJD(2);
	t320 = -t324 * qJD(1) + t325 * qJD(2);
	t319 = -t323 * qJD(1) + t341 * qJD(2);
	t318 = -t321 * t327 + t329 * t348;
	t317 = -t319 * t327 + t329 * t347;
	t316 = (t327 * t336 * t349 + (t329 * t331 + t332 * t362) * qJD(4)) * t330 + ((-t327 * t337 * t331 + t342 * t335) * qJD(4) + (t329 * t352 - t358) * t349 + ((-t329 * t358 + t352) * t331 - t333 * t362) * qJD(2)) * t328;
	t315 = (t321 * t359 + t322 * t336 + t332 * t346) * t331 - t318 * t335 + (t340 * t335 + (-t323 * t327 - t329 * t360) * t331) * qJD(4) + (-t324 * t332 + t344 * t336) * t349;
	t314 = (t319 * t359 + t320 * t336 + t332 * t345) * t331 - t317 * t335 + (t339 * t335 + (-t325 * t327 + t329 * t361) * t331) * qJD(4) + (t332 * t341 + t343 * t336) * t349;
	t1 = [0, t347, t317, t320 * t332 + (-t319 * t329 - t345) * t336 + t339 * qJD(3), t314, t314; 0, t348, t318, t322 * t332 + (-t321 * t329 - t346) * t336 + t340 * qJD(3), t315, t315; 0, 0, t328 * qJD(2) * t333 * t327, qJD(3) * t327 * t330 * t332 + (t342 * qJD(3) + (t329 * t356 + t357) * qJD(2)) * t328, t316, t316;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end