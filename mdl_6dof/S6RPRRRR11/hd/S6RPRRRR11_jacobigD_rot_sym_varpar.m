% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRRR11_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR11_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:13:37
	% EndTime: 2019-10-10 09:13:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:13:37
	% EndTime: 2019-10-10 09:13:37
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:13:37
	% EndTime: 2019-10-10 09:13:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:13:37
	% EndTime: 2019-10-10 09:13:37
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (18->12), div. (0->0), fcn. (18->8), ass. (0->7)
	t127 = sin(pkin(6)) * cos(pkin(7));
	t126 = cos(pkin(6)) * cos(pkin(13));
	t125 = cos(qJ(1));
	t124 = sin(qJ(1));
	t119 = sin(pkin(7));
	t118 = sin(pkin(13));
	t1 = [0, 0, (-(t118 * t124 - t125 * t126) * t119 + t125 * t127) * qJD(1), 0, 0, 0; 0, 0, (-(-t118 * t125 - t124 * t126) * t119 + t124 * t127) * qJD(1), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:13:38
	% EndTime: 2019-10-10 09:13:38
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (24->18), mult. (92->42), div. (0->0), fcn. (96->10), ass. (0->26)
	t185 = sin(pkin(7));
	t186 = sin(pkin(6));
	t206 = t186 * t185;
	t188 = cos(pkin(7));
	t192 = cos(qJ(3));
	t205 = t188 * t192;
	t184 = sin(pkin(13));
	t191 = sin(qJ(1));
	t204 = t191 * t184;
	t187 = cos(pkin(13));
	t203 = t191 * t187;
	t193 = cos(qJ(1));
	t202 = t193 * t184;
	t201 = t193 * t187;
	t200 = t191 * t206;
	t199 = t193 * t206;
	t198 = qJD(1) * t186 * t188;
	t189 = cos(pkin(6));
	t197 = t189 * t201 - t204;
	t196 = -t189 * t203 - t202;
	t195 = t189 * t202 + t203;
	t194 = -t189 * t204 + t201;
	t190 = sin(qJ(3));
	t183 = t196 * qJD(1);
	t182 = t197 * qJD(1);
	t1 = [0, 0, t182 * t185 + t193 * t198, t182 * t205 + (t194 * t192 + (t196 * t188 + t200) * t190) * qJD(3) + (-t195 * t190 - t192 * t199) * qJD(1), 0, 0; 0, 0, -t183 * t185 + t191 * t198, -t183 * t205 + (t195 * t192 + (t197 * t188 - t199) * t190) * qJD(3) + (t194 * t190 - t192 * t200) * qJD(1), 0, 0; 0, 0, 0, (t185 * t189 * t190 + (t187 * t188 * t190 + t184 * t192) * t186) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:13:38
	% EndTime: 2019-10-10 09:13:38
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (71->39), mult. (247->85), div. (0->0), fcn. (270->12), ass. (0->42)
	t262 = sin(pkin(7));
	t266 = cos(pkin(6));
	t291 = t262 * t266;
	t263 = sin(pkin(6));
	t269 = sin(qJ(1));
	t290 = t263 * t269;
	t272 = cos(qJ(1));
	t289 = t263 * t272;
	t265 = cos(pkin(7));
	t268 = sin(qJ(3));
	t288 = t265 * t268;
	t261 = sin(pkin(13));
	t287 = t269 * t261;
	t264 = cos(pkin(13));
	t286 = t269 * t264;
	t285 = t272 * t261;
	t284 = t272 * t264;
	t283 = qJD(1) * t263;
	t267 = sin(qJ(4));
	t282 = qJD(3) * t267;
	t281 = t269 * t283;
	t280 = t272 * t283;
	t279 = t262 * t281;
	t278 = t262 * t280;
	t257 = t266 * t284 - t287;
	t277 = t257 * t265 - t262 * t289;
	t259 = -t266 * t286 - t285;
	t276 = t259 * t265 + t262 * t290;
	t258 = t266 * t285 + t286;
	t260 = -t266 * t287 + t284;
	t271 = cos(qJ(3));
	t275 = t258 * t271 + t277 * t268;
	t274 = t260 * t271 + t276 * t268;
	t273 = t268 * t291 + (t261 * t271 + t264 * t288) * t263;
	t270 = cos(qJ(4));
	t256 = t260 * qJD(1);
	t255 = t259 * qJD(1);
	t254 = t258 * qJD(1);
	t253 = t257 * qJD(1);
	t252 = -t255 * t262 + t265 * t281;
	t251 = t253 * t262 + t265 * t280;
	t1 = [0, 0, t251, -t254 * t268 + (t253 * t265 - t278) * t271 + t274 * qJD(3), (-t253 * t288 - t254 * t271 + t268 * t278) * t267 - t251 * t270 + (t274 * t270 + (-t259 * t262 + t265 * t290) * t267) * qJD(4) + (-t260 * t268 + t276 * t271) * t282, 0; 0, 0, t252, t256 * t268 + (-t255 * t265 - t279) * t271 + t275 * qJD(3), (t255 * t288 + t256 * t271 + t268 * t279) * t267 - t252 * t270 + (t275 * t270 + (-t257 * t262 - t265 * t289) * t267) * qJD(4) + (-t258 * t268 + t277 * t271) * t282, 0; 0, 0, 0, t273 * qJD(3), (t273 * t270 + (-t263 * t264 * t262 + t266 * t265) * t267) * qJD(4) + (t271 * t291 + (t264 * t265 * t271 - t261 * t268) * t263) * t282, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:13:38
	% EndTime: 2019-10-10 09:13:38
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (118->39), mult. (402->85), div. (0->0), fcn. (444->12), ass. (0->45)
	t286 = sin(pkin(7));
	t290 = cos(pkin(6));
	t315 = t286 * t290;
	t287 = sin(pkin(6));
	t293 = sin(qJ(1));
	t314 = t287 * t293;
	t296 = cos(qJ(1));
	t313 = t287 * t296;
	t289 = cos(pkin(7));
	t292 = sin(qJ(3));
	t312 = t289 * t292;
	t285 = sin(pkin(13));
	t311 = t293 * t285;
	t288 = cos(pkin(13));
	t310 = t293 * t288;
	t309 = t296 * t285;
	t308 = t296 * t288;
	t307 = qJD(1) * t287;
	t291 = sin(qJ(4));
	t306 = qJD(3) * t291;
	t305 = t293 * t307;
	t304 = t296 * t307;
	t303 = t286 * t305;
	t302 = t286 * t304;
	t281 = t290 * t308 - t311;
	t301 = t281 * t289 - t286 * t313;
	t283 = -t290 * t310 - t309;
	t300 = t283 * t289 + t286 * t314;
	t282 = t290 * t309 + t310;
	t284 = -t290 * t311 + t308;
	t295 = cos(qJ(3));
	t299 = t282 * t295 + t301 * t292;
	t298 = t284 * t295 + t300 * t292;
	t297 = t292 * t315 + (t285 * t295 + t288 * t312) * t287;
	t294 = cos(qJ(4));
	t280 = t284 * qJD(1);
	t279 = t283 * qJD(1);
	t278 = t282 * qJD(1);
	t277 = t281 * qJD(1);
	t276 = -t279 * t286 + t289 * t305;
	t275 = t277 * t286 + t289 * t304;
	t274 = (t297 * t294 + (-t287 * t288 * t286 + t290 * t289) * t291) * qJD(4) + (t295 * t315 + (t288 * t289 * t295 - t285 * t292) * t287) * t306;
	t273 = (t279 * t312 + t280 * t295 + t292 * t303) * t291 - t276 * t294 + (t299 * t294 + (-t281 * t286 - t289 * t313) * t291) * qJD(4) + (-t282 * t292 + t301 * t295) * t306;
	t272 = (-t277 * t312 - t278 * t295 + t292 * t302) * t291 - t275 * t294 + (t298 * t294 + (-t283 * t286 + t289 * t314) * t291) * qJD(4) + (-t284 * t292 + t300 * t295) * t306;
	t1 = [0, 0, t275, -t278 * t292 + (t277 * t289 - t302) * t295 + t298 * qJD(3), t272, t272; 0, 0, t276, t280 * t292 + (-t279 * t289 - t303) * t295 + t299 * qJD(3), t273, t273; 0, 0, 0, t297 * qJD(3), t274, t274;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end