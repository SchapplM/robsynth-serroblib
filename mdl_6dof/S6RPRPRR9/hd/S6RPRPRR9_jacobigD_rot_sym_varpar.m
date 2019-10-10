% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRPRR9_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR9_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (18->12), div. (0->0), fcn. (18->8), ass. (0->7)
	t127 = sin(pkin(6)) * cos(pkin(7));
	t126 = cos(pkin(6)) * cos(pkin(12));
	t125 = cos(qJ(1));
	t124 = sin(qJ(1));
	t119 = sin(pkin(7));
	t118 = sin(pkin(12));
	t1 = [0, 0, (-(t118 * t124 - t125 * t126) * t119 + t125 * t127) * qJD(1), 0, 0, 0; 0, 0, (-(-t118 * t125 - t124 * t126) * t119 + t124 * t127) * qJD(1), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (18->12), div. (0->0), fcn. (18->8), ass. (0->7)
	t147 = sin(pkin(6)) * cos(pkin(7));
	t146 = cos(pkin(6)) * cos(pkin(12));
	t145 = cos(qJ(1));
	t144 = sin(qJ(1));
	t139 = sin(pkin(7));
	t138 = sin(pkin(12));
	t1 = [0, 0, (-(t138 * t144 - t145 * t146) * t139 + t145 * t147) * qJD(1), 0, 0, 0; 0, 0, (-(-t138 * t145 - t144 * t146) * t139 + t144 * t147) * qJD(1), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (39->20), mult. (141->45), div. (0->0), fcn. (151->12), ass. (0->34)
	t212 = sin(pkin(6));
	t218 = sin(qJ(1));
	t233 = t212 * t218;
	t220 = cos(qJ(1));
	t232 = t212 * t220;
	t210 = sin(pkin(12));
	t231 = t218 * t210;
	t214 = cos(pkin(12));
	t230 = t218 * t214;
	t229 = t220 * t210;
	t228 = t220 * t214;
	t215 = cos(pkin(7));
	t227 = qJD(1) * t212 * t215;
	t209 = sin(pkin(13));
	t213 = cos(pkin(13));
	t217 = sin(qJ(3));
	t219 = cos(qJ(3));
	t208 = -t219 * t209 - t217 * t213;
	t226 = t209 * t217 - t213 * t219;
	t216 = cos(pkin(6));
	t225 = t216 * t228 - t231;
	t224 = -t216 * t230 - t229;
	t223 = t216 * t229 + t230;
	t222 = t216 * t231 - t228;
	t221 = qJD(3) * t208;
	t211 = sin(pkin(7));
	t207 = t226 * qJD(3);
	t206 = t224 * qJD(1);
	t205 = t225 * qJD(1);
	t204 = t226 * t215;
	t203 = t226 * t211;
	t202 = t215 * t221;
	t201 = t211 * t221;
	t1 = [0, 0, t205 * t211 + t220 * t227, 0, t222 * t207 - t205 * t204 - t224 * t202 - t201 * t233 + (t203 * t232 + t223 * t208) * qJD(1), 0; 0, 0, -t206 * t211 + t218 * t227, 0, -t223 * t207 + t206 * t204 - t225 * t202 + t201 * t232 + (t203 * t233 + t222 * t208) * qJD(1), 0; 0, 0, 0, 0, -t216 * t201 + (-t202 * t214 - t207 * t210) * t212, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:11
	% EndTime: 2019-10-10 01:00:11
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (110->49), mult. (369->102), div. (0->0), fcn. (413->14), ass. (0->46)
	t292 = sin(pkin(13));
	t296 = cos(pkin(13));
	t301 = sin(qJ(3));
	t304 = cos(qJ(3));
	t307 = t301 * t292 - t304 * t296;
	t318 = t307 * qJD(3);
	t295 = sin(pkin(6));
	t302 = sin(qJ(1));
	t317 = t295 * t302;
	t305 = cos(qJ(1));
	t316 = t295 * t305;
	t293 = sin(pkin(12));
	t315 = t302 * t293;
	t297 = cos(pkin(12));
	t314 = t302 * t297;
	t313 = t305 * t293;
	t312 = t305 * t297;
	t311 = qJD(1) * t302;
	t310 = qJD(1) * t305;
	t298 = cos(pkin(7));
	t309 = qJD(1) * t295 * t298;
	t308 = t304 * t292 + t301 * t296;
	t299 = cos(pkin(6));
	t284 = t299 * t312 - t315;
	t286 = -t299 * t314 - t313;
	t285 = t299 * t313 + t314;
	t287 = -t299 * t315 + t312;
	t289 = t308 * qJD(3);
	t303 = cos(qJ(5));
	t300 = sin(qJ(5));
	t294 = sin(pkin(7));
	t283 = t287 * qJD(1);
	t282 = t286 * qJD(1);
	t281 = t285 * qJD(1);
	t280 = t284 * qJD(1);
	t279 = t308 * t298;
	t278 = t307 * t298;
	t277 = t308 * t294;
	t276 = t307 * t294;
	t275 = t298 * t318;
	t274 = t298 * t289;
	t273 = t294 * t318;
	t272 = t294 * t289;
	t271 = -t282 * t294 + t302 * t309;
	t270 = t280 * t294 + t305 * t309;
	t1 = [0, 0, t270, 0, t286 * t274 - t280 * t278 - t281 * t308 - t287 * t318 + (t272 * t302 + t276 * t310) * t295, (-t286 * t275 - t280 * t279 + t281 * t307 - t287 * t289 + (-t273 * t302 + t277 * t310) * t295) * t300 - t270 * t303 + ((t277 * t317 + t286 * t279 - t287 * t307) * t303 + (-t286 * t294 + t298 * t317) * t300) * qJD(5); 0, 0, t271, 0, t284 * t274 + t282 * t278 + t283 * t308 - t285 * t318 + (-t272 * t305 + t276 * t311) * t295, (-t284 * t275 + t282 * t279 - t283 * t307 - t285 * t289 + (t273 * t305 + t277 * t311) * t295) * t300 - t271 * t303 + ((-t277 * t316 + t284 * t279 - t285 * t307) * t303 + (-t284 * t294 - t298 * t316) * t300) * qJD(5); 0, 0, 0, 0, t299 * t272 + (t274 * t297 - t293 * t318) * t295, (-t299 * t273 + (-t275 * t297 - t289 * t293) * t295) * t300 + ((t299 * t277 + (t279 * t297 - t293 * t307) * t295) * t303 + (-t295 * t297 * t294 + t299 * t298) * t300) * qJD(5);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end