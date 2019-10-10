% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:21
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR5_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR5_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:28
	% EndTime: 2019-10-09 23:21:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:28
	% EndTime: 2019-10-09 23:21:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:28
	% EndTime: 2019-10-09 23:21:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:29
	% EndTime: 2019-10-09 23:21:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (15->10), div. (0->0), fcn. (15->7), ass. (0->7)
	t106 = sin(qJ(2));
	t109 = cos(pkin(6)) * t106;
	t108 = qJD(2) * sin(pkin(7));
	t107 = cos(qJ(2));
	t104 = cos(pkin(13));
	t102 = sin(pkin(13));
	t1 = [0, 0, -(t102 * t109 - t104 * t107) * t108, 0, 0, 0; 0, 0, -(-t102 * t107 - t104 * t109) * t108, 0, 0, 0; 0, 0, sin(pkin(6)) * t106 * t108, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:29
	% EndTime: 2019-10-09 23:21:29
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (24->17), mult. (88->43), div. (0->0), fcn. (92->10), ass. (0->24)
	t182 = sin(pkin(7));
	t183 = sin(pkin(6));
	t201 = t183 * t182;
	t185 = cos(pkin(7));
	t189 = cos(qJ(3));
	t200 = t185 * t189;
	t186 = cos(pkin(6));
	t188 = sin(qJ(2));
	t199 = t186 * t188;
	t190 = cos(qJ(2));
	t198 = t186 * t190;
	t187 = sin(qJ(3));
	t197 = t187 * t190;
	t196 = t188 * t189;
	t195 = qJD(2) * t187;
	t181 = sin(pkin(13));
	t184 = cos(pkin(13));
	t194 = -t181 * t188 + t184 * t198;
	t193 = t181 * t190 + t184 * t199;
	t192 = -t181 * t198 - t184 * t188;
	t191 = t181 * t199 - t184 * t190;
	t180 = t191 * qJD(2);
	t179 = t193 * qJD(2);
	t1 = [0, 0, -t180 * t182, -t180 * t200 + t192 * t195 + (-t191 * t189 + (t181 * t201 + t185 * t192) * t187) * qJD(3), 0, 0; 0, 0, t179 * t182, t179 * t200 + t194 * t195 + (t193 * t189 + (-t184 * t201 + t185 * t194) * t187) * qJD(3), 0, 0; 0, 0, qJD(2) * t188 * t201, t186 * t182 * qJD(3) * t187 + ((t185 * t197 + t196) * qJD(3) + (t185 * t196 + t197) * qJD(2)) * t183, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:29
	% EndTime: 2019-10-09 23:21:30
	% DurationCPUTime: 0.21s
	% Computational Cost: add. (70->40), mult. (240->91), div. (0->0), fcn. (263->12), ass. (0->38)
	t257 = sin(pkin(7));
	t258 = sin(pkin(6));
	t285 = t257 * t258;
	t265 = cos(qJ(4));
	t284 = t257 * t265;
	t260 = cos(pkin(7));
	t283 = t258 * t260;
	t263 = sin(qJ(3));
	t282 = t260 * t263;
	t266 = cos(qJ(3));
	t281 = t260 * t266;
	t261 = cos(pkin(6));
	t264 = sin(qJ(2));
	t280 = t261 * t264;
	t267 = cos(qJ(2));
	t279 = t261 * t267;
	t278 = t263 * t264;
	t277 = t263 * t267;
	t276 = t264 * t266;
	t275 = t266 * t267;
	t262 = sin(qJ(4));
	t274 = qJD(3) * t262;
	t256 = sin(pkin(13));
	t259 = cos(pkin(13));
	t252 = -t256 * t264 + t259 * t279;
	t273 = t252 * t260 - t259 * t285;
	t254 = -t256 * t279 - t259 * t264;
	t272 = t254 * t260 + t256 * t285;
	t253 = t256 * t267 + t259 * t280;
	t271 = t256 * t280 - t259 * t267;
	t270 = t260 * t277 + t276;
	t269 = t253 * t266 + t273 * t263;
	t268 = t272 * t263 - t266 * t271;
	t251 = t271 * qJD(2);
	t250 = t254 * qJD(2);
	t249 = t253 * qJD(2);
	t248 = t252 * qJD(2);
	t1 = [0, 0, -t251 * t257, t268 * qJD(3) + t250 * t263 - t251 * t281, (t250 * t266 + t251 * t282) * t262 + t251 * t284 + (t268 * t265 + (-t254 * t257 + t256 * t283) * t262) * qJD(4) + (t263 * t271 + t272 * t266) * t274, 0; 0, 0, t249 * t257, t269 * qJD(3) + t248 * t263 + t249 * t281, (t248 * t266 - t249 * t282) * t262 - t249 * t284 + (t269 * t265 + (-t252 * t257 - t259 * t283) * t262) * qJD(4) + (-t253 * t263 + t273 * t266) * t274, 0; 0, 0, qJD(2) * t264 * t285, t261 * t257 * qJD(3) * t263 + (t270 * qJD(3) + (t260 * t276 + t277) * qJD(2)) * t258, (t257 * t266 * t274 + (t260 * t262 + t263 * t284) * qJD(4)) * t261 + ((-t257 * t267 * t262 + t270 * t265) * qJD(4) + (t260 * t275 - t278) * t274 + ((-t260 * t278 + t275) * t262 - t264 * t284) * qJD(2)) * t258, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:21:29
	% EndTime: 2019-10-09 23:21:30
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (116->40), mult. (392->91), div. (0->0), fcn. (434->12), ass. (0->41)
	t285 = sin(pkin(7));
	t286 = sin(pkin(6));
	t313 = t285 * t286;
	t293 = cos(qJ(4));
	t312 = t285 * t293;
	t288 = cos(pkin(7));
	t311 = t286 * t288;
	t291 = sin(qJ(3));
	t310 = t288 * t291;
	t294 = cos(qJ(3));
	t309 = t288 * t294;
	t289 = cos(pkin(6));
	t292 = sin(qJ(2));
	t308 = t289 * t292;
	t295 = cos(qJ(2));
	t307 = t289 * t295;
	t306 = t291 * t292;
	t305 = t291 * t295;
	t304 = t292 * t294;
	t303 = t294 * t295;
	t290 = sin(qJ(4));
	t302 = qJD(3) * t290;
	t284 = sin(pkin(13));
	t287 = cos(pkin(13));
	t280 = -t284 * t292 + t287 * t307;
	t301 = t280 * t288 - t287 * t313;
	t282 = -t284 * t307 - t287 * t292;
	t300 = t282 * t288 + t284 * t313;
	t281 = t284 * t295 + t287 * t308;
	t299 = t284 * t308 - t287 * t295;
	t298 = t288 * t305 + t304;
	t297 = t281 * t294 + t301 * t291;
	t296 = t300 * t291 - t294 * t299;
	t279 = t299 * qJD(2);
	t278 = t282 * qJD(2);
	t277 = t281 * qJD(2);
	t276 = t280 * qJD(2);
	t275 = (t285 * t294 * t302 + (t288 * t290 + t291 * t312) * qJD(4)) * t289 + ((-t285 * t295 * t290 + t298 * t293) * qJD(4) + (t288 * t303 - t306) * t302 + ((-t288 * t306 + t303) * t290 - t292 * t312) * qJD(2)) * t286;
	t274 = (t278 * t294 + t279 * t310) * t290 + t279 * t312 + (t296 * t293 + (-t282 * t285 + t284 * t311) * t290) * qJD(4) + (t291 * t299 + t300 * t294) * t302;
	t273 = (t276 * t294 - t277 * t310) * t290 - t277 * t312 + (t297 * t293 + (-t280 * t285 - t287 * t311) * t290) * qJD(4) + (-t281 * t291 + t301 * t294) * t302;
	t1 = [0, 0, -t279 * t285, t296 * qJD(3) + t278 * t291 - t279 * t309, t274, t274; 0, 0, t277 * t285, t297 * qJD(3) + t276 * t291 + t277 * t309, t273, t273; 0, 0, qJD(2) * t292 * t313, t289 * t285 * qJD(3) * t291 + (t298 * qJD(3) + (t288 * t304 + t305) * qJD(2)) * t286, t275, t275;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end