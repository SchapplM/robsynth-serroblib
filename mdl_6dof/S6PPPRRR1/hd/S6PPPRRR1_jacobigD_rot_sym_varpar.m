% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPPRRR1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPPRRR1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobigD_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:22
	% EndTime: 2019-10-10 08:49:22
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (31->23), mult. (101->55), div. (0->0), fcn. (126->14), ass. (0->27)
	t208 = sin(pkin(12));
	t217 = cos(pkin(6));
	t227 = t208 * t217;
	t210 = sin(pkin(7));
	t211 = sin(pkin(6));
	t226 = t210 * t211;
	t225 = t210 * t217;
	t216 = cos(pkin(7));
	t224 = t211 * t216;
	t213 = cos(pkin(13));
	t223 = t213 * t216;
	t214 = cos(pkin(12));
	t222 = t214 * t217;
	t207 = sin(pkin(13));
	t202 = -t208 * t207 + t213 * t222;
	t221 = t202 * t216 - t214 * t226;
	t204 = -t214 * t207 - t213 * t227;
	t220 = t204 * t216 + t208 * t226;
	t219 = cos(qJ(4));
	t218 = sin(qJ(4));
	t215 = cos(pkin(8));
	t212 = cos(pkin(14));
	t209 = sin(pkin(8));
	t206 = sin(pkin(14));
	t205 = -t207 * t227 + t214 * t213;
	t203 = t207 * t222 + t208 * t213;
	t1 = [0, 0, 0, 0, ((t205 * t212 + t220 * t206) * t219 + ((-t205 * t206 + t220 * t212) * t215 + (-t204 * t210 + t208 * t224) * t209) * t218) * qJD(4), 0; 0, 0, 0, 0, ((t203 * t212 + t221 * t206) * t219 + ((-t203 * t206 + t221 * t212) * t215 + (-t202 * t210 - t214 * t224) * t209) * t218) * qJD(4), 0; 0, 0, 0, 0, ((t211 * t207 * t212 + (t211 * t223 + t225) * t206) * t219 + ((t212 * t225 + (-t206 * t207 + t212 * t223) * t211) * t215 + (-t213 * t226 + t217 * t216) * t209) * t218) * qJD(4), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:23
	% EndTime: 2019-10-10 08:49:23
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (117->35), mult. (361->80), div. (0->0), fcn. (460->16), ass. (0->45)
	t297 = sin(pkin(12));
	t306 = cos(pkin(6));
	t325 = t297 * t306;
	t299 = sin(pkin(7));
	t300 = sin(pkin(6));
	t324 = t299 * t300;
	t323 = t299 * t306;
	t305 = cos(pkin(7));
	t322 = t300 * t305;
	t302 = cos(pkin(13));
	t321 = t302 * t305;
	t303 = cos(pkin(12));
	t320 = t303 * t306;
	t307 = sin(qJ(5));
	t319 = qJD(4) * t307;
	t296 = sin(pkin(13));
	t292 = t296 * t320 + t297 * t302;
	t295 = sin(pkin(14));
	t301 = cos(pkin(14));
	t291 = -t297 * t296 + t302 * t320;
	t315 = t291 * t305 - t303 * t324;
	t282 = -t292 * t295 + t315 * t301;
	t288 = -t291 * t299 - t303 * t322;
	t298 = sin(pkin(8));
	t304 = cos(pkin(8));
	t318 = t282 * t304 + t288 * t298;
	t294 = -t296 * t325 + t303 * t302;
	t293 = -t303 * t296 - t302 * t325;
	t314 = t293 * t305 + t297 * t324;
	t284 = -t294 * t295 + t314 * t301;
	t289 = -t293 * t299 + t297 * t322;
	t317 = t284 * t304 + t289 * t298;
	t286 = t301 * t323 + (-t295 * t296 + t301 * t321) * t300;
	t290 = -t302 * t324 + t306 * t305;
	t316 = t286 * t304 + t290 * t298;
	t283 = t292 * t301 + t315 * t295;
	t308 = sin(qJ(4));
	t310 = cos(qJ(4));
	t313 = t283 * t310 + t318 * t308;
	t285 = t294 * t301 + t314 * t295;
	t312 = t285 * t310 + t317 * t308;
	t287 = t300 * t296 * t301 + (t300 * t321 + t323) * t295;
	t311 = t287 * t310 + t316 * t308;
	t309 = cos(qJ(5));
	t1 = [0, 0, 0, 0, t312 * qJD(4), (t312 * t309 + (-t284 * t298 + t289 * t304) * t307) * qJD(5) + (-t285 * t308 + t317 * t310) * t319; 0, 0, 0, 0, t313 * qJD(4), (t313 * t309 + (-t282 * t298 + t288 * t304) * t307) * qJD(5) + (-t283 * t308 + t318 * t310) * t319; 0, 0, 0, 0, t311 * qJD(4), (t311 * t309 + (-t286 * t298 + t290 * t304) * t307) * qJD(5) + (-t287 * t308 + t316 * t310) * t319;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end