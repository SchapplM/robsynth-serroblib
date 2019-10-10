% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR11
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
% Datum: 2019-10-10 01:04
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRPRR11_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR11_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
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
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:00
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (4->4), mult. (18->12), div. (0->0), fcn. (18->8), ass. (0->7)
	t168 = sin(pkin(6)) * cos(pkin(7));
	t167 = cos(pkin(6)) * cos(pkin(12));
	t166 = cos(qJ(1));
	t165 = sin(qJ(1));
	t160 = sin(pkin(7));
	t159 = sin(pkin(12));
	t1 = [0, 0, (-(t159 * t165 - t166 * t167) * t160 + t166 * t168) * qJD(1), 0, 0, 0; 0, 0, (-(-t159 * t166 - t165 * t167) * t160 + t165 * t168) * qJD(1), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:00
	% EndTime: 2019-10-10 01:04:01
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (24->18), mult. (92->42), div. (0->0), fcn. (96->10), ass. (0->26)
	t193 = sin(pkin(7));
	t194 = sin(pkin(6));
	t214 = t194 * t193;
	t196 = cos(pkin(7));
	t200 = cos(qJ(3));
	t213 = t196 * t200;
	t192 = sin(pkin(12));
	t199 = sin(qJ(1));
	t212 = t199 * t192;
	t195 = cos(pkin(12));
	t211 = t199 * t195;
	t201 = cos(qJ(1));
	t210 = t201 * t192;
	t209 = t201 * t195;
	t208 = t199 * t214;
	t207 = t201 * t214;
	t206 = qJD(1) * t194 * t196;
	t197 = cos(pkin(6));
	t205 = t197 * t209 - t212;
	t204 = -t197 * t211 - t210;
	t203 = t197 * t210 + t211;
	t202 = -t197 * t212 + t209;
	t198 = sin(qJ(3));
	t191 = t204 * qJD(1);
	t190 = t205 * qJD(1);
	t1 = [0, 0, t190 * t193 + t201 * t206, 0, t190 * t213 + (t202 * t200 + (t204 * t196 + t208) * t198) * qJD(3) + (-t203 * t198 - t200 * t207) * qJD(1), 0; 0, 0, -t191 * t193 + t199 * t206, 0, -t191 * t213 + (t203 * t200 + (t205 * t196 - t207) * t198) * qJD(3) + (t202 * t198 - t200 * t208) * qJD(1), 0; 0, 0, 0, 0, (t193 * t197 * t198 + (t195 * t196 * t198 + t192 * t200) * t194) * qJD(3), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:04:01
	% EndTime: 2019-10-10 01:04:01
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (82->40), mult. (247->85), div. (0->0), fcn. (270->12), ass. (0->43)
	t275 = sin(pkin(7));
	t279 = cos(pkin(6));
	t302 = t275 * t279;
	t276 = sin(pkin(6));
	t281 = sin(qJ(1));
	t301 = t276 * t281;
	t283 = cos(qJ(1));
	t300 = t276 * t283;
	t278 = cos(pkin(7));
	t280 = sin(qJ(3));
	t299 = t278 * t280;
	t274 = sin(pkin(12));
	t298 = t281 * t274;
	t277 = cos(pkin(12));
	t297 = t281 * t277;
	t296 = t283 * t274;
	t295 = t283 * t277;
	t294 = qJD(1) * t276;
	t273 = pkin(13) + qJ(5);
	t271 = sin(t273);
	t293 = qJD(3) * t271;
	t292 = t281 * t294;
	t291 = t283 * t294;
	t290 = t275 * t292;
	t289 = t275 * t291;
	t267 = t279 * t295 - t298;
	t288 = t267 * t278 - t275 * t300;
	t269 = -t279 * t297 - t296;
	t287 = t269 * t278 + t275 * t301;
	t268 = t279 * t296 + t297;
	t270 = -t279 * t298 + t295;
	t282 = cos(qJ(3));
	t286 = t268 * t282 + t288 * t280;
	t285 = t270 * t282 + t287 * t280;
	t284 = t280 * t302 + (t274 * t282 + t277 * t299) * t276;
	t272 = cos(t273);
	t266 = t270 * qJD(1);
	t265 = t269 * qJD(1);
	t264 = t268 * qJD(1);
	t263 = t267 * qJD(1);
	t262 = -t265 * t275 + t278 * t292;
	t261 = t263 * t275 + t278 * t291;
	t1 = [0, 0, t261, 0, -t264 * t280 + (t263 * t278 - t289) * t282 + t285 * qJD(3), (-t263 * t299 - t264 * t282 + t280 * t289) * t271 - t261 * t272 + (t285 * t272 + (-t269 * t275 + t278 * t301) * t271) * qJD(5) + (-t270 * t280 + t287 * t282) * t293; 0, 0, t262, 0, t266 * t280 + (-t265 * t278 - t290) * t282 + t286 * qJD(3), (t265 * t299 + t266 * t282 + t280 * t290) * t271 - t262 * t272 + (t286 * t272 + (-t267 * t275 - t278 * t300) * t271) * qJD(5) + (-t268 * t280 + t288 * t282) * t293; 0, 0, 0, 0, t284 * qJD(3), (t284 * t272 + (-t276 * t277 * t275 + t279 * t278) * t271) * qJD(5) + (t282 * t302 + (t277 * t278 * t282 - t274 * t280) * t276) * t293;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end