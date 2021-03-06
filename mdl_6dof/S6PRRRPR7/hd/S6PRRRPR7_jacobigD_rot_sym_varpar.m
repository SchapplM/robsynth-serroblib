% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR7_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR7_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.03s
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
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:05
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
	t181 = sin(pkin(12));
	t184 = cos(pkin(12));
	t194 = -t181 * t188 + t184 * t198;
	t193 = t181 * t190 + t184 * t199;
	t192 = -t181 * t198 - t184 * t188;
	t191 = t181 * t199 - t184 * t190;
	t180 = t191 * qJD(2);
	t179 = t193 * qJD(2);
	t1 = [0, 0, -t180 * t182, -t180 * t200 + t192 * t195 + (-t191 * t189 + (t181 * t201 + t192 * t185) * t187) * qJD(3), 0, 0; 0, 0, t179 * t182, t179 * t200 + t194 * t195 + (t193 * t189 + (-t184 * t201 + t194 * t185) * t187) * qJD(3), 0, 0; 0, 0, qJD(2) * t188 * t201, t186 * t182 * qJD(3) * t187 + ((t185 * t197 + t196) * qJD(3) + (t185 * t196 + t197) * qJD(2)) * t183, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:05
	% EndTime: 2019-10-09 22:58:05
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (24->17), mult. (88->43), div. (0->0), fcn. (92->10), ass. (0->24)
	t230 = sin(pkin(7));
	t231 = sin(pkin(6));
	t249 = t231 * t230;
	t233 = cos(pkin(7));
	t237 = cos(qJ(3));
	t248 = t233 * t237;
	t234 = cos(pkin(6));
	t236 = sin(qJ(2));
	t247 = t234 * t236;
	t238 = cos(qJ(2));
	t246 = t234 * t238;
	t235 = sin(qJ(3));
	t245 = t235 * t238;
	t244 = t236 * t237;
	t243 = qJD(2) * t235;
	t229 = sin(pkin(12));
	t232 = cos(pkin(12));
	t242 = -t229 * t236 + t232 * t246;
	t241 = t229 * t238 + t232 * t247;
	t240 = -t229 * t246 - t232 * t236;
	t239 = t229 * t247 - t232 * t238;
	t228 = t239 * qJD(2);
	t227 = t241 * qJD(2);
	t1 = [0, 0, -t228 * t230, -t228 * t248 + t240 * t243 + (-t239 * t237 + (t229 * t249 + t240 * t233) * t235) * qJD(3), 0, 0; 0, 0, t227 * t230, t227 * t248 + t242 * t243 + (t241 * t237 + (-t232 * t249 + t242 * t233) * t235) * qJD(3), 0, 0; 0, 0, qJD(2) * t236 * t249, t234 * t230 * qJD(3) * t235 + ((t233 * t245 + t244) * qJD(3) + (t233 * t244 + t245) * qJD(2)) * t231, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:05
	% EndTime: 2019-10-09 22:58:05
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (70->40), mult. (240->91), div. (0->0), fcn. (263->12), ass. (0->38)
	t266 = sin(pkin(7));
	t267 = sin(pkin(6));
	t294 = t266 * t267;
	t274 = cos(qJ(4));
	t293 = t266 * t274;
	t269 = cos(pkin(7));
	t292 = t267 * t269;
	t272 = sin(qJ(3));
	t291 = t269 * t272;
	t275 = cos(qJ(3));
	t290 = t269 * t275;
	t270 = cos(pkin(6));
	t273 = sin(qJ(2));
	t289 = t270 * t273;
	t276 = cos(qJ(2));
	t288 = t270 * t276;
	t287 = t272 * t273;
	t286 = t272 * t276;
	t285 = t273 * t275;
	t284 = t275 * t276;
	t271 = sin(qJ(4));
	t283 = qJD(3) * t271;
	t265 = sin(pkin(12));
	t268 = cos(pkin(12));
	t261 = -t265 * t273 + t268 * t288;
	t282 = t261 * t269 - t268 * t294;
	t263 = -t265 * t288 - t268 * t273;
	t281 = t263 * t269 + t265 * t294;
	t262 = t265 * t276 + t268 * t289;
	t280 = t265 * t289 - t268 * t276;
	t279 = t269 * t286 + t285;
	t278 = t262 * t275 + t282 * t272;
	t277 = t281 * t272 - t275 * t280;
	t260 = t280 * qJD(2);
	t259 = t263 * qJD(2);
	t258 = t262 * qJD(2);
	t257 = t261 * qJD(2);
	t1 = [0, 0, -t260 * t266, t277 * qJD(3) + t259 * t272 - t260 * t290, 0, (t259 * t275 + t260 * t291) * t271 + t260 * t293 + (t277 * t274 + (-t263 * t266 + t265 * t292) * t271) * qJD(4) + (t272 * t280 + t281 * t275) * t283; 0, 0, t258 * t266, t278 * qJD(3) + t257 * t272 + t258 * t290, 0, (t257 * t275 - t258 * t291) * t271 - t258 * t293 + (t278 * t274 + (-t261 * t266 - t268 * t292) * t271) * qJD(4) + (-t262 * t272 + t282 * t275) * t283; 0, 0, qJD(2) * t273 * t294, t270 * t266 * qJD(3) * t272 + (t279 * qJD(3) + (t269 * t285 + t286) * qJD(2)) * t267, 0, (t266 * t275 * t283 + (t269 * t271 + t272 * t293) * qJD(4)) * t270 + ((-t266 * t276 * t271 + t279 * t274) * qJD(4) + (t269 * t284 - t287) * t283 + ((-t269 * t287 + t284) * t271 - t273 * t293) * qJD(2)) * t267;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end