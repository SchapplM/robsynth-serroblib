% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:11
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPRRRR10_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRR10_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:29
	% EndTime: 2019-10-10 09:11:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (18->12), div. (0->0), fcn. (18->8), ass. (0->7)
	t127 = sin(pkin(6)) * cos(pkin(7));
	t126 = cos(pkin(6)) * cos(pkin(13));
	t125 = cos(qJ(1));
	t124 = sin(qJ(1));
	t119 = sin(pkin(7));
	t118 = sin(pkin(13));
	t1 = [0, 0, (-(t118 * t124 - t125 * t126) * t119 + t125 * t127) * qJD(1), 0, 0, 0; 0, 0, (-(-t118 * t125 - t124 * t126) * t119 + t124 * t127) * qJD(1), 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (24->18), mult. (92->42), div. (0->0), fcn. (96->10), ass. (0->26)
	t185 = sin(pkin(7));
	t186 = sin(pkin(6));
	t206 = t186 * t185;
	t188 = cos(pkin(7));
	t192 = cos(qJ(3));
	t205 = t188 * t192;
	t184 = sin(pkin(13));
	t191 = sin(qJ(1));
	t204 = t191 * t184;
	t187 = cos(pkin(13));
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
	% StartTime: 2019-10-10 09:11:30
	% EndTime: 2019-10-10 09:11:30
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (44->18), mult. (166->42), div. (0->0), fcn. (174->10), ass. (0->29)
	t210 = sin(pkin(7));
	t211 = sin(pkin(6));
	t231 = t211 * t210;
	t213 = cos(pkin(7));
	t217 = cos(qJ(3));
	t230 = t213 * t217;
	t209 = sin(pkin(13));
	t216 = sin(qJ(1));
	t229 = t216 * t209;
	t212 = cos(pkin(13));
	t228 = t216 * t212;
	t218 = cos(qJ(1));
	t227 = t218 * t209;
	t226 = t218 * t212;
	t225 = t216 * t231;
	t224 = t218 * t231;
	t223 = qJD(1) * t211 * t213;
	t214 = cos(pkin(6));
	t222 = t214 * t226 - t229;
	t221 = -t214 * t228 - t227;
	t220 = t214 * t227 + t228;
	t219 = -t214 * t229 + t226;
	t215 = sin(qJ(3));
	t208 = t221 * qJD(1);
	t207 = t222 * qJD(1);
	t206 = (t210 * t214 * t215 + (t212 * t213 * t215 + t209 * t217) * t211) * qJD(3);
	t205 = -t208 * t230 + (t220 * t217 + (t222 * t213 - t224) * t215) * qJD(3) + (t219 * t215 - t217 * t225) * qJD(1);
	t204 = t207 * t230 + (t219 * t217 + (t221 * t213 + t225) * t215) * qJD(3) + (-t220 * t215 - t217 * t224) * qJD(1);
	t1 = [0, 0, t207 * t210 + t218 * t223, t204, t204, 0; 0, 0, -t208 * t210 + t216 * t223, t205, t205, 0; 0, 0, 0, t206, t206, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:11:31
	% EndTime: 2019-10-10 09:11:31
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (108->39), mult. (321->80), div. (0->0), fcn. (348->12), ass. (0->47)
	t303 = sin(pkin(7));
	t307 = cos(pkin(6));
	t330 = t303 * t307;
	t304 = sin(pkin(6));
	t309 = sin(qJ(1));
	t329 = t304 * t309;
	t311 = cos(qJ(1));
	t328 = t304 * t311;
	t305 = cos(pkin(13));
	t306 = cos(pkin(7));
	t327 = t305 * t306;
	t302 = sin(pkin(13));
	t326 = t309 * t302;
	t325 = t309 * t305;
	t324 = t311 * t302;
	t323 = t311 * t305;
	t322 = qJD(1) * t304;
	t301 = qJ(4) + qJ(5);
	t298 = sin(t301);
	t321 = qJD(3) * t298;
	t320 = t309 * t322;
	t319 = t311 * t322;
	t294 = t307 * t323 - t326;
	t318 = t294 * t306 - t303 * t328;
	t296 = -t307 * t325 - t324;
	t317 = t296 * t306 + t303 * t329;
	t295 = t307 * t324 + t325;
	t297 = -t307 * t326 + t323;
	t290 = t294 * qJD(1);
	t316 = -t290 * t306 + t303 * t319;
	t292 = t296 * qJD(1);
	t315 = t292 * t306 + t303 * t320;
	t308 = sin(qJ(3));
	t310 = cos(qJ(3));
	t314 = t295 * t310 + t318 * t308;
	t313 = t297 * t310 + t317 * t308;
	t312 = t308 * t330 + (t302 * t310 + t308 * t327) * t304;
	t300 = qJD(4) + qJD(5);
	t299 = cos(t301);
	t293 = t297 * qJD(1);
	t291 = t295 * qJD(1);
	t289 = -t292 * t303 + t306 * t320;
	t288 = t290 * t303 + t306 * t319;
	t287 = t312 * qJD(3);
	t286 = t314 * qJD(3) + t293 * t308 - t315 * t310;
	t285 = t313 * qJD(3) - t291 * t308 - t316 * t310;
	t1 = [0, 0, t288, t285, t285, (t313 * t300 - t288) * t299 + (-t291 * t310 + (-t296 * t303 + t306 * t329) * t300 + t316 * t308) * t298 + (-t297 * t308 + t317 * t310) * t321; 0, 0, t289, t286, t286, (t314 * t300 - t289) * t299 + (t293 * t310 + (-t294 * t303 - t306 * t328) * t300 + t315 * t308) * t298 + (-t295 * t308 + t318 * t310) * t321; 0, 0, 0, t287, t287, (t312 * t299 + (-t304 * t305 * t303 + t307 * t306) * t298) * t300 + (t310 * t330 + (-t302 * t308 + t310 * t327) * t304) * t321;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end