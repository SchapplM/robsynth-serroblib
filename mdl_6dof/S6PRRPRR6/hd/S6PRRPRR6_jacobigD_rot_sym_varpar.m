% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR6
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
% Datum: 2019-10-09 22:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR6_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR6_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:20
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
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:20
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (15->10), div. (0->0), fcn. (15->7), ass. (0->7)
	t170 = sin(qJ(2));
	t173 = cos(pkin(6)) * t170;
	t172 = qJD(2) * sin(pkin(7));
	t171 = cos(qJ(2));
	t168 = cos(pkin(12));
	t166 = sin(pkin(12));
	t1 = [0, 0, -(t166 * t173 - t168 * t171) * t172, 0, 0, 0; 0, 0, -(-t166 * t171 - t168 * t173) * t172, 0, 0, 0; 0, 0, sin(pkin(6)) * t170 * t172, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:20
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (24->17), mult. (88->43), div. (0->0), fcn. (92->10), ass. (0->24)
	t190 = sin(pkin(7));
	t191 = sin(pkin(6));
	t209 = t191 * t190;
	t193 = cos(pkin(7));
	t197 = cos(qJ(3));
	t208 = t193 * t197;
	t194 = cos(pkin(6));
	t196 = sin(qJ(2));
	t207 = t194 * t196;
	t198 = cos(qJ(2));
	t206 = t194 * t198;
	t195 = sin(qJ(3));
	t205 = t195 * t198;
	t204 = t196 * t197;
	t203 = qJD(2) * t195;
	t189 = sin(pkin(12));
	t192 = cos(pkin(12));
	t202 = -t189 * t196 + t192 * t206;
	t201 = t189 * t198 + t192 * t207;
	t200 = -t189 * t206 - t192 * t196;
	t199 = t189 * t207 - t192 * t198;
	t188 = t199 * qJD(2);
	t187 = t201 * qJD(2);
	t1 = [0, 0, -t188 * t190, 0, -t188 * t208 + t200 * t203 + (-t199 * t197 + (t189 * t209 + t200 * t193) * t195) * qJD(3), 0; 0, 0, t187 * t190, 0, t187 * t208 + t202 * t203 + (t201 * t197 + (-t192 * t209 + t202 * t193) * t195) * qJD(3), 0; 0, 0, qJD(2) * t196 * t209, 0, t194 * t190 * qJD(3) * t195 + ((t193 * t205 + t204) * qJD(3) + (t193 * t204 + t205) * qJD(2)) * t191, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:21
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (82->41), mult. (240->92), div. (0->0), fcn. (263->12), ass. (0->41)
	t266 = pkin(13) + qJ(5);
	t265 = cos(t266);
	t268 = sin(pkin(7));
	t296 = t268 * t265;
	t269 = sin(pkin(6));
	t295 = t268 * t269;
	t274 = sin(qJ(2));
	t294 = t268 * t274;
	t271 = cos(pkin(7));
	t293 = t269 * t271;
	t273 = sin(qJ(3));
	t292 = t271 * t273;
	t275 = cos(qJ(3));
	t291 = t271 * t275;
	t272 = cos(pkin(6));
	t290 = t272 * t274;
	t276 = cos(qJ(2));
	t289 = t272 * t276;
	t288 = t273 * t274;
	t287 = t273 * t276;
	t286 = t274 * t275;
	t285 = t275 * t276;
	t264 = sin(t266);
	t284 = qJD(3) * t264;
	t283 = qJD(3) * t268;
	t267 = sin(pkin(12));
	t270 = cos(pkin(12));
	t260 = -t267 * t274 + t270 * t289;
	t282 = t260 * t271 - t270 * t295;
	t262 = -t267 * t289 - t270 * t274;
	t281 = t262 * t271 + t267 * t295;
	t261 = t267 * t276 + t270 * t290;
	t280 = t267 * t290 - t270 * t276;
	t279 = t271 * t287 + t286;
	t278 = t261 * t275 + t282 * t273;
	t277 = t281 * t273 - t275 * t280;
	t259 = t280 * qJD(2);
	t258 = t262 * qJD(2);
	t257 = t261 * qJD(2);
	t256 = t260 * qJD(2);
	t1 = [0, 0, -t259 * t268, 0, t277 * qJD(3) + t258 * t273 - t259 * t291, (t258 * t275 + t259 * t292) * t264 + t259 * t296 + (t277 * t265 + (-t262 * t268 + t267 * t293) * t264) * qJD(5) + (t273 * t280 + t281 * t275) * t284; 0, 0, t257 * t268, 0, t278 * qJD(3) + t256 * t273 + t257 * t291, (t256 * t275 - t257 * t292) * t264 - t257 * t296 + (t278 * t265 + (-t260 * t268 - t270 * t293) * t264) * qJD(5) + (-t261 * t273 + t282 * t275) * t284; 0, 0, t269 * qJD(2) * t294, 0, t272 * t273 * t283 + (t279 * qJD(3) + (t271 * t286 + t287) * qJD(2)) * t269, (t275 * t264 * t283 + (t264 * t271 + t273 * t296) * qJD(5)) * t272 + ((-t268 * t276 * t264 + t279 * t265) * qJD(5) + (t271 * t285 - t288) * t284 + ((-t271 * t288 + t285) * t264 - t265 * t294) * qJD(2)) * t269;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end