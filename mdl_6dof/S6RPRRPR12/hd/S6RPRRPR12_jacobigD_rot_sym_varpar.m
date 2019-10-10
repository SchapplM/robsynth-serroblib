% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:43
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRPR12_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR12_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:26
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
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:26
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (24->18), mult. (92->42), div. (0->0), fcn. (96->10), ass. (0->26)
	t185 = sin(pkin(7));
	t186 = sin(pkin(6));
	t206 = t186 * t185;
	t188 = cos(pkin(7));
	t192 = cos(qJ(3));
	t205 = t188 * t192;
	t184 = sin(pkin(12));
	t191 = sin(qJ(1));
	t204 = t191 * t184;
	t187 = cos(pkin(12));
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
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:26
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (24->18), mult. (92->42), div. (0->0), fcn. (96->10), ass. (0->26)
	t205 = sin(pkin(7));
	t206 = sin(pkin(6));
	t226 = t206 * t205;
	t208 = cos(pkin(7));
	t212 = cos(qJ(3));
	t225 = t208 * t212;
	t204 = sin(pkin(12));
	t211 = sin(qJ(1));
	t224 = t211 * t204;
	t207 = cos(pkin(12));
	t223 = t211 * t207;
	t213 = cos(qJ(1));
	t222 = t213 * t204;
	t221 = t213 * t207;
	t220 = t211 * t226;
	t219 = t213 * t226;
	t218 = qJD(1) * t206 * t208;
	t209 = cos(pkin(6));
	t217 = t209 * t221 - t224;
	t216 = -t209 * t223 - t222;
	t215 = t209 * t222 + t223;
	t214 = -t209 * t224 + t221;
	t210 = sin(qJ(3));
	t203 = t216 * qJD(1);
	t202 = t217 * qJD(1);
	t1 = [0, 0, t202 * t205 + t213 * t218, t202 * t225 + (t214 * t212 + (t216 * t208 + t220) * t210) * qJD(3) + (-t215 * t210 - t212 * t219) * qJD(1), 0, 0; 0, 0, -t203 * t205 + t211 * t218, -t203 * t225 + (t215 * t212 + (t217 * t208 - t219) * t210) * qJD(3) + (t214 * t210 - t212 * t220) * qJD(1), 0, 0; 0, 0, 0, (t205 * t209 * t210 + (t207 * t208 * t210 + t204 * t212) * t206) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:27
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (71->39), mult. (247->85), div. (0->0), fcn. (270->12), ass. (0->42)
	t270 = sin(pkin(7));
	t274 = cos(pkin(6));
	t299 = t270 * t274;
	t271 = sin(pkin(6));
	t277 = sin(qJ(1));
	t298 = t271 * t277;
	t280 = cos(qJ(1));
	t297 = t271 * t280;
	t273 = cos(pkin(7));
	t276 = sin(qJ(3));
	t296 = t273 * t276;
	t269 = sin(pkin(12));
	t295 = t277 * t269;
	t272 = cos(pkin(12));
	t294 = t277 * t272;
	t293 = t280 * t269;
	t292 = t280 * t272;
	t291 = qJD(1) * t271;
	t278 = cos(qJ(4));
	t290 = qJD(3) * t278;
	t289 = t277 * t291;
	t288 = t280 * t291;
	t287 = t270 * t289;
	t286 = t270 * t288;
	t265 = t274 * t292 - t295;
	t285 = t265 * t273 - t270 * t297;
	t267 = -t274 * t294 - t293;
	t284 = t267 * t273 + t270 * t298;
	t266 = t274 * t293 + t294;
	t268 = -t274 * t295 + t292;
	t279 = cos(qJ(3));
	t283 = t266 * t279 + t285 * t276;
	t282 = t268 * t279 + t284 * t276;
	t281 = t276 * t299 + (t269 * t279 + t272 * t296) * t271;
	t275 = sin(qJ(4));
	t264 = t268 * qJD(1);
	t263 = t267 * qJD(1);
	t262 = t266 * qJD(1);
	t261 = t265 * qJD(1);
	t260 = -t263 * t270 + t273 * t289;
	t259 = t261 * t270 + t273 * t288;
	t1 = [0, 0, t259, -t262 * t276 + (t261 * t273 - t286) * t279 + t282 * qJD(3), 0, (-t261 * t296 - t262 * t279 + t276 * t286) * t278 + t259 * t275 + (-t282 * t275 + (-t267 * t270 + t273 * t298) * t278) * qJD(4) + (-t268 * t276 + t284 * t279) * t290; 0, 0, t260, t264 * t276 + (-t263 * t273 - t287) * t279 + t283 * qJD(3), 0, (t263 * t296 + t264 * t279 + t276 * t287) * t278 + t260 * t275 + (-t283 * t275 + (-t265 * t270 - t273 * t297) * t278) * qJD(4) + (-t266 * t276 + t285 * t279) * t290; 0, 0, 0, t281 * qJD(3), 0, (-t281 * t275 + (-t271 * t272 * t270 + t274 * t273) * t278) * qJD(4) + (t279 * t299 + (t272 * t273 * t279 - t269 * t276) * t271) * t290;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end