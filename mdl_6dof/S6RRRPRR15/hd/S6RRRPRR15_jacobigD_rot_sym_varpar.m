% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR15_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR15_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (8->8), mult. (35->23), div. (0->0), fcn. (35->8), ass. (0->15)
	t139 = sin(pkin(6));
	t152 = t139 * cos(pkin(7));
	t142 = sin(qJ(2));
	t143 = sin(qJ(1));
	t151 = t142 * t143;
	t145 = cos(qJ(1));
	t150 = t142 * t145;
	t144 = cos(qJ(2));
	t149 = t143 * t144;
	t148 = t144 * t145;
	t147 = qJD(1) * t139;
	t138 = sin(pkin(7));
	t146 = qJD(2) * t138;
	t141 = cos(pkin(6));
	t1 = [0, t145 * t147, -(t141 * t151 - t148) * t146 + (-(-t141 * t148 + t151) * t138 + t145 * t152) * qJD(1), 0, 0, 0; 0, t143 * t147, -(-t141 * t150 - t149) * t146 + (-(-t141 * t149 - t150) * t138 + t143 * t152) * qJD(1), 0, 0, 0; 0, 0, t139 * t142 * t146, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (8->8), mult. (35->23), div. (0->0), fcn. (35->8), ass. (0->15)
	t166 = sin(pkin(6));
	t179 = t166 * cos(pkin(7));
	t169 = sin(qJ(2));
	t170 = sin(qJ(1));
	t178 = t169 * t170;
	t172 = cos(qJ(1));
	t177 = t169 * t172;
	t171 = cos(qJ(2));
	t176 = t170 * t171;
	t175 = t171 * t172;
	t174 = qJD(1) * t166;
	t165 = sin(pkin(7));
	t173 = qJD(2) * t165;
	t168 = cos(pkin(6));
	t1 = [0, t172 * t174, -(t168 * t178 - t175) * t173 + (-(-t168 * t175 + t178) * t165 + t172 * t179) * qJD(1), 0, 0, 0; 0, t170 * t174, -(-t168 * t177 - t176) * t173 + (-(-t168 * t176 - t177) * t165 + t170 * t179) * qJD(1), 0, 0, 0; 0, 0, t166 * t169 * t173, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:24
	% EndTime: 2019-10-10 12:18:24
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (38->24), mult. (136->54), div. (0->0), fcn. (140->10), ass. (0->31)
	t215 = sin(pkin(7));
	t216 = sin(pkin(6));
	t242 = t216 * t215;
	t217 = cos(pkin(7));
	t219 = sin(qJ(3));
	t241 = t217 * t219;
	t220 = sin(qJ(2));
	t240 = t219 * t220;
	t221 = sin(qJ(1));
	t239 = t221 * t220;
	t223 = cos(qJ(2));
	t238 = t221 * t223;
	t222 = cos(qJ(3));
	t237 = t222 * t223;
	t224 = cos(qJ(1));
	t236 = t224 * t220;
	t235 = t224 * t223;
	t234 = qJD(1) * t216;
	t233 = qJD(2) * t222;
	t232 = t221 * t242;
	t231 = t224 * t242;
	t230 = t221 * t234;
	t229 = t224 * t234;
	t218 = cos(pkin(6));
	t228 = t218 * t235 - t239;
	t227 = -t218 * t238 - t236;
	t226 = t218 * t236 + t238;
	t225 = t218 * t239 - t235;
	t214 = t227 * qJD(1) - t226 * qJD(2);
	t213 = -t228 * qJD(1) + t225 * qJD(2);
	t1 = [0, t229, -t213 * t215 + t217 * t229, 0, t213 * t241 + t227 * t233 + (t219 * t231 - t226 * t222) * qJD(1) + (t225 * t219 + (t227 * t217 + t232) * t222) * qJD(3), 0; 0, t230, -t214 * t215 + t217 * t230, 0, t214 * t241 + t228 * t233 + (t219 * t232 - t225 * t222) * qJD(1) + (-t226 * t219 + (t228 * t217 - t231) * t222) * qJD(3), 0; 0, 0, qJD(2) * t220 * t242, 0, t218 * t215 * qJD(3) * t222 + ((t217 * t237 - t240) * qJD(3) + (-t217 * t240 + t237) * qJD(2)) * t216, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:25
	% EndTime: 2019-10-10 12:18:25
	% DurationCPUTime: 0.34s
	% Computational Cost: add. (100->49), mult. (332->103), div. (0->0), fcn. (355->12), ass. (0->47)
	t307 = sin(qJ(3));
	t311 = cos(qJ(3));
	t305 = cos(pkin(6));
	t312 = cos(qJ(2));
	t313 = cos(qJ(1));
	t324 = t313 * t312;
	t308 = sin(qJ(2));
	t309 = sin(qJ(1));
	t328 = t309 * t308;
	t314 = t305 * t328 - t324;
	t325 = t313 * t308;
	t327 = t309 * t312;
	t300 = -t305 * t327 - t325;
	t302 = sin(pkin(7));
	t304 = cos(pkin(7));
	t303 = sin(pkin(6));
	t334 = t303 * t309;
	t316 = t300 * t304 + t302 * t334;
	t340 = t307 * t314 + t316 * t311;
	t299 = t305 * t325 + t327;
	t298 = t305 * t324 - t328;
	t333 = t303 * t313;
	t317 = -t298 * t304 + t302 * t333;
	t339 = t299 * t307 + t317 * t311;
	t336 = t302 * t308;
	t335 = t302 * t311;
	t332 = t304 * t311;
	t331 = t307 * t308;
	t330 = t307 * t312;
	t329 = t308 * t311;
	t326 = t311 * t312;
	t323 = qJD(1) * t303;
	t310 = cos(qJ(5));
	t322 = qJD(3) * t310;
	t321 = t309 * t323;
	t320 = t313 * t323;
	t319 = t302 * t321;
	t318 = t302 * t320;
	t315 = t304 * t326 - t331;
	t306 = sin(qJ(5));
	t297 = -t314 * qJD(1) + t298 * qJD(2);
	t296 = t300 * qJD(1) - t299 * qJD(2);
	t295 = -t299 * qJD(1) + t300 * qJD(2);
	t294 = -t298 * qJD(1) + t314 * qJD(2);
	t293 = -t296 * t302 + t304 * t321;
	t292 = -t294 * t302 + t304 * t320;
	t1 = [0, t320, t292, 0, t295 * t311 + (t294 * t304 + t318) * t307 + t340 * qJD(3), t292 * t306 - (-t294 * t332 + t295 * t307 - t311 * t318) * t310 + ((-t300 * t302 + t304 * t334) * t310 - t340 * t306) * qJD(5) - (t316 * t307 - t311 * t314) * t322; 0, t321, t293, 0, t297 * t311 + (t296 * t304 + t319) * t307 - t339 * qJD(3), t293 * t306 - (-t296 * t332 + t297 * t307 - t311 * t319) * t310 + ((-t298 * t302 - t304 * t333) * t310 + t339 * t306) * qJD(5) - (t299 * t311 - t317 * t307) * t322; 0, 0, t303 * qJD(2) * t336, 0, t305 * qJD(3) * t335 + (t315 * qJD(3) + (-t304 * t331 + t326) * qJD(2)) * t303, (-t302 * t307 * t322 + (t304 * t310 - t306 * t335) * qJD(5)) * t305 + ((-t302 * t312 * t310 - t315 * t306) * qJD(5) - (t304 * t330 + t329) * t322 + (t306 * t336 - (t304 * t329 + t330) * t310) * qJD(2)) * t303;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end