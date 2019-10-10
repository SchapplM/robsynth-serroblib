% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR5
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
% Datum: 2019-10-09 22:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR5_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR5_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:08
	% EndTime: 2019-10-09 22:54:08
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
	% StartTime: 2019-10-09 22:54:08
	% EndTime: 2019-10-09 22:54:08
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
	% StartTime: 2019-10-09 22:54:08
	% EndTime: 2019-10-09 22:54:08
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (24->17), mult. (88->43), div. (0->0), fcn. (92->10), ass. (0->24)
	t193 = sin(pkin(7));
	t194 = sin(pkin(6));
	t212 = t194 * t193;
	t196 = cos(pkin(7));
	t200 = cos(qJ(3));
	t211 = t196 * t200;
	t197 = cos(pkin(6));
	t199 = sin(qJ(2));
	t210 = t197 * t199;
	t201 = cos(qJ(2));
	t209 = t197 * t201;
	t198 = sin(qJ(3));
	t208 = t198 * t201;
	t207 = t199 * t200;
	t206 = qJD(2) * t198;
	t192 = sin(pkin(12));
	t195 = cos(pkin(12));
	t205 = -t192 * t199 + t195 * t209;
	t204 = t192 * t201 + t195 * t210;
	t203 = -t192 * t209 - t195 * t199;
	t202 = t192 * t210 - t195 * t201;
	t191 = t202 * qJD(2);
	t190 = t204 * qJD(2);
	t1 = [0, 0, -t191 * t193, -t191 * t211 + t203 * t206 + (-t202 * t200 + (t192 * t212 + t203 * t196) * t198) * qJD(3), 0, 0; 0, 0, t190 * t193, t190 * t211 + t205 * t206 + (t204 * t200 + (-t195 * t212 + t205 * t196) * t198) * qJD(3), 0, 0; 0, 0, qJD(2) * t199 * t212, t197 * t193 * qJD(3) * t198 + ((t196 * t208 + t207) * qJD(3) + (t196 * t207 + t208) * qJD(2)) * t194, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:09
	% EndTime: 2019-10-09 22:54:09
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (82->41), mult. (240->92), div. (0->0), fcn. (263->12), ass. (0->41)
	t267 = qJ(4) + pkin(13);
	t266 = cos(t267);
	t269 = sin(pkin(7));
	t297 = t269 * t266;
	t270 = sin(pkin(6));
	t296 = t269 * t270;
	t275 = sin(qJ(2));
	t295 = t269 * t275;
	t272 = cos(pkin(7));
	t294 = t270 * t272;
	t274 = sin(qJ(3));
	t293 = t272 * t274;
	t276 = cos(qJ(3));
	t292 = t272 * t276;
	t273 = cos(pkin(6));
	t291 = t273 * t275;
	t277 = cos(qJ(2));
	t290 = t273 * t277;
	t289 = t274 * t275;
	t288 = t274 * t277;
	t287 = t275 * t276;
	t286 = t276 * t277;
	t265 = sin(t267);
	t285 = qJD(3) * t265;
	t284 = qJD(3) * t269;
	t268 = sin(pkin(12));
	t271 = cos(pkin(12));
	t261 = -t268 * t275 + t271 * t290;
	t283 = t261 * t272 - t271 * t296;
	t263 = -t268 * t290 - t271 * t275;
	t282 = t263 * t272 + t268 * t296;
	t262 = t268 * t277 + t271 * t291;
	t281 = t268 * t291 - t271 * t277;
	t280 = t272 * t288 + t287;
	t279 = t262 * t276 + t283 * t274;
	t278 = t282 * t274 - t276 * t281;
	t260 = t281 * qJD(2);
	t259 = t263 * qJD(2);
	t258 = t262 * qJD(2);
	t257 = t261 * qJD(2);
	t1 = [0, 0, -t260 * t269, t278 * qJD(3) + t259 * t274 - t260 * t292, 0, (t259 * t276 + t260 * t293) * t265 + t260 * t297 + (t278 * t266 + (-t263 * t269 + t268 * t294) * t265) * qJD(4) + (t274 * t281 + t282 * t276) * t285; 0, 0, t258 * t269, t279 * qJD(3) + t257 * t274 + t258 * t292, 0, (t257 * t276 - t258 * t293) * t265 - t258 * t297 + (t279 * t266 + (-t261 * t269 - t271 * t294) * t265) * qJD(4) + (-t262 * t274 + t283 * t276) * t285; 0, 0, t270 * qJD(2) * t295, t273 * t274 * t284 + (t280 * qJD(3) + (t272 * t287 + t288) * qJD(2)) * t270, 0, (t276 * t265 * t284 + (t265 * t272 + t274 * t297) * qJD(4)) * t273 + ((-t269 * t277 * t265 + t280 * t266) * qJD(4) + (t272 * t286 - t289) * t285 + ((-t272 * t289 + t286) * t265 - t266 * t295) * qJD(2)) * t270;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end