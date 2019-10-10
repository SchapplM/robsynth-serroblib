% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:19
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR4_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR4_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (15->10), div. (0->0), fcn. (15->7), ass. (0->7)
	t106 = sin(qJ(2));
	t109 = cos(pkin(6)) * t106;
	t108 = qJD(2) * sin(pkin(7));
	t107 = cos(qJ(2));
	t104 = cos(pkin(13));
	t102 = sin(pkin(13));
	t1 = [0, 0, -(t102 * t109 - t104 * t107) * t108, 0, 0, 0; 0, 0, -(-t102 * t107 - t104 * t109) * t108, 0, 0, 0; 0, 0, sin(pkin(6)) * t106 * t108, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:24
	% EndTime: 2019-10-09 23:19:24
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
	t181 = sin(pkin(13));
	t184 = cos(pkin(13));
	t194 = -t181 * t188 + t184 * t198;
	t193 = t181 * t190 + t184 * t199;
	t192 = -t181 * t198 - t184 * t188;
	t191 = t181 * t199 - t184 * t190;
	t180 = t191 * qJD(2);
	t179 = t193 * qJD(2);
	t1 = [0, 0, -t180 * t182, -t180 * t200 + t192 * t195 + (-t191 * t189 + (t181 * t201 + t185 * t192) * t187) * qJD(3), 0, 0; 0, 0, t179 * t182, t179 * t200 + t194 * t195 + (t193 * t189 + (-t184 * t201 + t185 * t194) * t187) * qJD(3), 0, 0; 0, 0, qJD(2) * t188 * t201, t186 * t182 * qJD(3) * t187 + ((t185 * t197 + t196) * qJD(3) + (t185 * t196 + t197) * qJD(2)) * t183, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:24
	% EndTime: 2019-10-09 23:19:24
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (44->17), mult. (161->43), div. (0->0), fcn. (169->10), ass. (0->27)
	t209 = sin(pkin(7));
	t210 = sin(pkin(6));
	t228 = t210 * t209;
	t212 = cos(pkin(7));
	t216 = cos(qJ(3));
	t227 = t212 * t216;
	t213 = cos(pkin(6));
	t215 = sin(qJ(2));
	t226 = t213 * t215;
	t217 = cos(qJ(2));
	t225 = t213 * t217;
	t214 = sin(qJ(3));
	t224 = t214 * t217;
	t223 = t215 * t216;
	t222 = qJD(2) * t214;
	t208 = sin(pkin(13));
	t211 = cos(pkin(13));
	t221 = -t208 * t215 + t211 * t225;
	t220 = t208 * t217 + t211 * t226;
	t219 = -t208 * t225 - t211 * t215;
	t218 = t208 * t226 - t211 * t217;
	t207 = t218 * qJD(2);
	t206 = t220 * qJD(2);
	t205 = t213 * t209 * qJD(3) * t214 + ((t212 * t224 + t223) * qJD(3) + (t212 * t223 + t224) * qJD(2)) * t210;
	t204 = -t207 * t227 + t219 * t222 + (-t218 * t216 + (t208 * t228 + t219 * t212) * t214) * qJD(3);
	t203 = t206 * t227 + t221 * t222 + (t220 * t216 + (-t211 * t228 + t221 * t212) * t214) * qJD(3);
	t1 = [0, 0, -t207 * t209, t204, t204, 0; 0, 0, t206 * t209, t203, t203, 0; 0, 0, qJD(2) * t215 * t228, t205, t205, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:25
	% EndTime: 2019-10-09 23:19:25
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (108->42), mult. (313->90), div. (0->0), fcn. (340->12), ass. (0->46)
	t296 = sin(pkin(13));
	t299 = cos(pkin(13));
	t305 = cos(qJ(2));
	t301 = cos(pkin(6));
	t303 = sin(qJ(2));
	t319 = t301 * t303;
	t289 = t296 * t305 + t299 * t319;
	t285 = t289 * qJD(2);
	t297 = sin(pkin(7));
	t326 = t285 * t297;
	t309 = t296 * t319 - t299 * t305;
	t287 = t309 * qJD(2);
	t325 = t287 * t297;
	t298 = sin(pkin(6));
	t324 = t297 * t298;
	t323 = t297 * t303;
	t300 = cos(pkin(7));
	t322 = t298 * t300;
	t302 = sin(qJ(3));
	t321 = t300 * t302;
	t304 = cos(qJ(3));
	t320 = t300 * t304;
	t318 = t301 * t305;
	t317 = t302 * t303;
	t316 = t302 * t305;
	t315 = t303 * t304;
	t314 = t304 * t305;
	t295 = qJ(4) + qJ(5);
	t292 = sin(t295);
	t313 = qJD(3) * t292;
	t312 = qJD(3) * t297;
	t288 = -t296 * t303 + t299 * t318;
	t311 = t288 * t300 - t299 * t324;
	t290 = -t296 * t318 - t299 * t303;
	t310 = t290 * t300 + t296 * t324;
	t308 = t300 * t316 + t315;
	t307 = t289 * t304 + t311 * t302;
	t306 = t310 * t302 - t304 * t309;
	t294 = qJD(4) + qJD(5);
	t293 = cos(t295);
	t286 = t290 * qJD(2);
	t284 = t288 * qJD(2);
	t283 = t301 * t302 * t312 + (t308 * qJD(3) + (t300 * t315 + t316) * qJD(2)) * t298;
	t282 = t306 * qJD(3) + t286 * t302 - t287 * t320;
	t281 = t307 * qJD(3) + t284 * t302 + t285 * t320;
	t1 = [0, 0, -t325, t282, t282, (t306 * t294 + t325) * t293 + (t287 * t321 + t286 * t304 + (-t290 * t297 + t296 * t322) * t294) * t292 + (t302 * t309 + t310 * t304) * t313; 0, 0, t326, t281, t281, (t307 * t294 - t326) * t293 + (-t285 * t321 + t284 * t304 + (-t288 * t297 - t299 * t322) * t294) * t292 + (-t289 * t302 + t311 * t304) * t313; 0, 0, t298 * qJD(2) * t323, t283, t283, (t297 * t302 * t294 * t293 + (t300 * t294 + t304 * t312) * t292) * t301 + ((-t297 * t305 * t292 + t308 * t293) * t294 + (t300 * t314 - t317) * t313 + ((-t300 * t317 + t314) * t292 - t293 * t323) * qJD(2)) * t298;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end