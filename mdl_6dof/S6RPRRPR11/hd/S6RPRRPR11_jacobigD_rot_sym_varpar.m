% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR11
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
% Datum: 2019-10-10 01:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRPR11_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR11_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
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
	% StartTime: 2019-10-10 01:41:21
	% EndTime: 2019-10-10 01:41:21
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
	% StartTime: 2019-10-10 01:41:21
	% EndTime: 2019-10-10 01:41:21
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (24->18), mult. (92->42), div. (0->0), fcn. (96->10), ass. (0->26)
	t225 = sin(pkin(7));
	t226 = sin(pkin(6));
	t246 = t226 * t225;
	t228 = cos(pkin(7));
	t232 = cos(qJ(3));
	t245 = t228 * t232;
	t224 = sin(pkin(12));
	t231 = sin(qJ(1));
	t244 = t231 * t224;
	t227 = cos(pkin(12));
	t243 = t231 * t227;
	t233 = cos(qJ(1));
	t242 = t233 * t224;
	t241 = t233 * t227;
	t240 = t231 * t246;
	t239 = t233 * t246;
	t238 = qJD(1) * t226 * t228;
	t229 = cos(pkin(6));
	t237 = t229 * t241 - t244;
	t236 = -t229 * t243 - t242;
	t235 = t229 * t242 + t243;
	t234 = -t229 * t244 + t241;
	t230 = sin(qJ(3));
	t223 = t236 * qJD(1);
	t222 = t237 * qJD(1);
	t1 = [0, 0, t222 * t225 + t233 * t238, t222 * t245 + (t234 * t232 + (t236 * t228 + t240) * t230) * qJD(3) + (-t235 * t230 - t232 * t239) * qJD(1), 0, 0; 0, 0, -t223 * t225 + t231 * t238, -t223 * t245 + (t235 * t232 + (t237 * t228 - t239) * t230) * qJD(3) + (t234 * t230 - t232 * t240) * qJD(1), 0, 0; 0, 0, 0, (t225 * t229 * t230 + (t227 * t228 * t230 + t224 * t232) * t226) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:21
	% EndTime: 2019-10-10 01:41:22
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (71->39), mult. (247->85), div. (0->0), fcn. (270->12), ass. (0->42)
	t275 = sin(pkin(7));
	t279 = cos(pkin(6));
	t304 = t275 * t279;
	t276 = sin(pkin(6));
	t282 = sin(qJ(1));
	t303 = t276 * t282;
	t285 = cos(qJ(1));
	t302 = t276 * t285;
	t278 = cos(pkin(7));
	t281 = sin(qJ(3));
	t301 = t278 * t281;
	t274 = sin(pkin(12));
	t300 = t282 * t274;
	t277 = cos(pkin(12));
	t299 = t282 * t277;
	t298 = t285 * t274;
	t297 = t285 * t277;
	t296 = qJD(1) * t276;
	t280 = sin(qJ(4));
	t295 = qJD(3) * t280;
	t294 = t282 * t296;
	t293 = t285 * t296;
	t292 = t275 * t294;
	t291 = t275 * t293;
	t270 = t279 * t297 - t300;
	t290 = t270 * t278 - t275 * t302;
	t272 = -t279 * t299 - t298;
	t289 = t272 * t278 + t275 * t303;
	t271 = t279 * t298 + t299;
	t273 = -t279 * t300 + t297;
	t284 = cos(qJ(3));
	t288 = t271 * t284 + t290 * t281;
	t287 = t273 * t284 + t289 * t281;
	t286 = t281 * t304 + (t274 * t284 + t277 * t301) * t276;
	t283 = cos(qJ(4));
	t269 = t273 * qJD(1);
	t268 = t272 * qJD(1);
	t267 = t271 * qJD(1);
	t266 = t270 * qJD(1);
	t265 = -t268 * t275 + t278 * t294;
	t264 = t266 * t275 + t278 * t293;
	t1 = [0, 0, t264, -t267 * t281 + (t266 * t278 - t291) * t284 + t287 * qJD(3), 0, (-t266 * t301 - t267 * t284 + t281 * t291) * t280 - t264 * t283 + (t287 * t283 + (-t272 * t275 + t278 * t303) * t280) * qJD(4) + (-t273 * t281 + t289 * t284) * t295; 0, 0, t265, t269 * t281 + (-t268 * t278 - t292) * t284 + t288 * qJD(3), 0, (t268 * t301 + t269 * t284 + t281 * t292) * t280 - t265 * t283 + (t288 * t283 + (-t270 * t275 - t278 * t302) * t280) * qJD(4) + (-t271 * t281 + t290 * t284) * t295; 0, 0, 0, t286 * qJD(3), 0, (t286 * t283 + (-t276 * t277 * t275 + t279 * t278) * t280) * qJD(4) + (t284 * t304 + (t277 * t278 * t284 - t274 * t281) * t276) * t295;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end