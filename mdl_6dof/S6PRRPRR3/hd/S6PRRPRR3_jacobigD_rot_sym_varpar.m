% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR3_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR3_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (15->10), div. (0->0), fcn. (15->7), ass. (0->7)
	t106 = sin(qJ(2));
	t109 = cos(pkin(6)) * t106;
	t108 = qJD(2) * sin(pkin(7));
	t107 = cos(qJ(2));
	t104 = cos(pkin(12));
	t102 = sin(pkin(12));
	t1 = [0, 0, -(t102 * t109 - t104 * t107) * t108, 0, 0, 0; 0, 0, -(-t102 * t107 - t104 * t109) * t108, 0, 0, 0; 0, 0, sin(pkin(6)) * t106 * t108, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (15->10), div. (0->0), fcn. (15->7), ass. (0->7)
	t129 = sin(qJ(2));
	t132 = cos(pkin(6)) * t129;
	t131 = qJD(2) * sin(pkin(7));
	t130 = cos(qJ(2));
	t127 = cos(pkin(12));
	t125 = sin(pkin(12));
	t1 = [0, 0, -(t125 * t132 - t127 * t130) * t131, 0, 0, 0; 0, 0, -(-t125 * t130 - t127 * t132) * t131, 0, 0, 0; 0, 0, sin(pkin(6)) * t129 * t131, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:27
	% EndTime: 2019-10-09 22:29:27
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (39->18), mult. (136->43), div. (0->0), fcn. (146->12), ass. (0->28)
	t210 = cos(pkin(7));
	t204 = sin(pkin(13));
	t208 = cos(pkin(13));
	t212 = sin(qJ(3));
	t214 = cos(qJ(3));
	t203 = -t214 * t204 - t212 * t208;
	t216 = qJD(3) * t203;
	t198 = t210 * t216;
	t226 = qJD(2) * t203 + t198;
	t206 = sin(pkin(7));
	t197 = t206 * t216;
	t207 = sin(pkin(6));
	t225 = t207 * t197;
	t211 = cos(pkin(6));
	t213 = sin(qJ(2));
	t224 = t211 * t213;
	t215 = cos(qJ(2));
	t223 = t211 * t215;
	t221 = t204 * t212 - t208 * t214;
	t205 = sin(pkin(12));
	t209 = cos(pkin(12));
	t219 = t205 * t215 + t209 * t224;
	t217 = t205 * t224 - t209 * t215;
	t202 = t221 * qJD(3);
	t201 = t217 * qJD(2);
	t200 = t219 * qJD(2);
	t199 = t221 * t210;
	t1 = [0, 0, -t201 * t206, 0, t201 * t199 + t217 * t202 - t205 * t225 + t226 * (t205 * t223 + t209 * t213), 0; 0, 0, t200 * t206, 0, -t200 * t199 - t219 * t202 + t209 * t225 - t226 * (-t205 * t213 + t209 * t223), 0; 0, 0, t207 * qJD(2) * t213 * t206, 0, -t211 * t197 + (-t198 * t215 - t202 * t213 + (-t199 * t213 - t203 * t215) * qJD(2)) * t207, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:27
	% EndTime: 2019-10-09 22:29:28
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (109->49), mult. (360->103), div. (0->0), fcn. (404->14), ass. (0->40)
	t288 = sin(pkin(13));
	t292 = cos(pkin(13));
	t297 = sin(qJ(3));
	t300 = cos(qJ(3));
	t304 = t297 * t288 - t300 * t292;
	t312 = t304 * qJD(3);
	t289 = sin(pkin(12));
	t291 = sin(pkin(6));
	t311 = t289 * t291;
	t290 = sin(pkin(7));
	t299 = cos(qJ(5));
	t310 = t290 * t299;
	t293 = cos(pkin(12));
	t309 = t291 * t293;
	t294 = cos(pkin(7));
	t308 = t291 * t294;
	t295 = cos(pkin(6));
	t298 = sin(qJ(2));
	t307 = t295 * t298;
	t301 = cos(qJ(2));
	t306 = t295 * t301;
	t305 = t300 * t288 + t297 * t292;
	t280 = -t289 * t298 + t293 * t306;
	t281 = t289 * t301 + t293 * t307;
	t282 = -t289 * t306 - t293 * t298;
	t303 = t289 * t307 - t293 * t301;
	t285 = t305 * qJD(3);
	t296 = sin(qJ(5));
	t279 = t303 * qJD(2);
	t278 = t282 * qJD(2);
	t277 = t281 * qJD(2);
	t276 = t280 * qJD(2);
	t275 = t305 * t294;
	t274 = t304 * t294;
	t273 = t305 * t290;
	t272 = t294 * t312;
	t271 = t294 * t285;
	t270 = t290 * t312;
	t269 = t290 * t285;
	t1 = [0, 0, -t279 * t290, 0, t269 * t311 + t282 * t271 + t279 * t274 + t278 * t305 + t303 * t312, (-t270 * t311 - t282 * t272 + t279 * t275 - t278 * t304 + t285 * t303) * t296 + t279 * t310 + ((t273 * t311 + t282 * t275 + t303 * t304) * t299 + (-t282 * t290 + t289 * t308) * t296) * qJD(5); 0, 0, t277 * t290, 0, -t269 * t309 + t280 * t271 - t277 * t274 + t276 * t305 - t281 * t312, (t270 * t309 - t280 * t272 - t277 * t275 - t276 * t304 - t281 * t285) * t296 - t277 * t310 + ((-t273 * t309 + t280 * t275 - t281 * t304) * t299 + (-t280 * t290 - t293 * t308) * t296) * qJD(5); 0, 0, t291 * qJD(2) * t298 * t290, 0, t295 * t269 + (t271 * t301 - t312 * t298 + (-t274 * t298 + t301 * t305) * qJD(2)) * t291, (-t270 * t296 + (t273 * t299 + t294 * t296) * qJD(5)) * t295 + ((-t272 * t301 - t285 * t298) * t296 + ((t275 * t301 - t298 * t304) * t299 - t290 * t301 * t296) * qJD(5) + ((-t275 * t298 - t301 * t304) * t296 - t298 * t310) * qJD(2)) * t291;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end