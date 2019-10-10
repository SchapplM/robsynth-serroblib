% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:37
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRPR9_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR9_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:31
	% EndTime: 2019-10-10 01:37:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
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
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
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
	% StartTime: 2019-10-10 01:37:32
	% EndTime: 2019-10-10 01:37:32
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (24->18), mult. (92->42), div. (0->0), fcn. (96->10), ass. (0->26)
	t196 = sin(pkin(7));
	t197 = sin(pkin(6));
	t217 = t197 * t196;
	t199 = cos(pkin(7));
	t203 = cos(qJ(3));
	t216 = t199 * t203;
	t195 = sin(pkin(12));
	t202 = sin(qJ(1));
	t215 = t202 * t195;
	t198 = cos(pkin(12));
	t214 = t202 * t198;
	t204 = cos(qJ(1));
	t213 = t204 * t195;
	t212 = t204 * t198;
	t211 = t202 * t217;
	t210 = t204 * t217;
	t209 = qJD(1) * t197 * t199;
	t200 = cos(pkin(6));
	t208 = t200 * t212 - t215;
	t207 = -t200 * t214 - t213;
	t206 = t200 * t213 + t214;
	t205 = -t200 * t215 + t212;
	t201 = sin(qJ(3));
	t194 = t207 * qJD(1);
	t193 = t208 * qJD(1);
	t1 = [0, 0, t193 * t196 + t204 * t209, t193 * t216 + (t205 * t203 + (t207 * t199 + t211) * t201) * qJD(3) + (-t206 * t201 - t203 * t210) * qJD(1), 0, 0; 0, 0, -t194 * t196 + t202 * t209, -t194 * t216 + (t206 * t203 + (t208 * t199 - t210) * t201) * qJD(3) + (t205 * t201 - t203 * t211) * qJD(1), 0, 0; 0, 0, 0, (t196 * t200 * t201 + (t198 * t199 * t201 + t195 * t203) * t197) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:37:33
	% EndTime: 2019-10-10 01:37:33
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (82->40), mult. (247->85), div. (0->0), fcn. (270->12), ass. (0->43)
	t274 = sin(pkin(7));
	t278 = cos(pkin(6));
	t301 = t274 * t278;
	t275 = sin(pkin(6));
	t280 = sin(qJ(1));
	t300 = t275 * t280;
	t282 = cos(qJ(1));
	t299 = t275 * t282;
	t277 = cos(pkin(7));
	t279 = sin(qJ(3));
	t298 = t277 * t279;
	t273 = sin(pkin(12));
	t297 = t280 * t273;
	t276 = cos(pkin(12));
	t296 = t280 * t276;
	t295 = t282 * t273;
	t294 = t282 * t276;
	t293 = qJD(1) * t275;
	t272 = qJ(4) + pkin(13);
	t270 = sin(t272);
	t292 = qJD(3) * t270;
	t291 = t280 * t293;
	t290 = t282 * t293;
	t289 = t274 * t291;
	t288 = t274 * t290;
	t266 = t278 * t294 - t297;
	t287 = t266 * t277 - t274 * t299;
	t268 = -t278 * t296 - t295;
	t286 = t268 * t277 + t274 * t300;
	t267 = t278 * t295 + t296;
	t269 = -t278 * t297 + t294;
	t281 = cos(qJ(3));
	t285 = t267 * t281 + t287 * t279;
	t284 = t269 * t281 + t286 * t279;
	t283 = t279 * t301 + (t273 * t281 + t276 * t298) * t275;
	t271 = cos(t272);
	t265 = t269 * qJD(1);
	t264 = t268 * qJD(1);
	t263 = t267 * qJD(1);
	t262 = t266 * qJD(1);
	t261 = -t264 * t274 + t277 * t291;
	t260 = t262 * t274 + t277 * t290;
	t1 = [0, 0, t260, -t263 * t279 + (t262 * t277 - t288) * t281 + t284 * qJD(3), 0, (-t262 * t298 - t263 * t281 + t279 * t288) * t270 - t260 * t271 + (t284 * t271 + (-t268 * t274 + t277 * t300) * t270) * qJD(4) + (-t269 * t279 + t286 * t281) * t292; 0, 0, t261, t265 * t279 + (-t264 * t277 - t289) * t281 + t285 * qJD(3), 0, (t264 * t298 + t265 * t281 + t279 * t289) * t270 - t261 * t271 + (t285 * t271 + (-t266 * t274 - t277 * t299) * t270) * qJD(4) + (-t267 * t279 + t287 * t281) * t292; 0, 0, 0, t283 * qJD(3), 0, (t283 * t271 + (-t275 * t276 * t274 + t278 * t277) * t270) * qJD(4) + (t281 * t301 + (t276 * t277 * t281 - t273 * t279) * t275) * t292;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end