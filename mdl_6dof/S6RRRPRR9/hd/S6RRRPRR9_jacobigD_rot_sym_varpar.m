% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:08
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR9_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR9_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:08
	% EndTime: 2019-10-10 12:08:08
	% DurationCPUTime: 0.06s
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
	% StartTime: 2019-10-10 12:08:08
	% EndTime: 2019-10-10 12:08:08
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (8->8), mult. (35->23), div. (0->0), fcn. (35->8), ass. (0->15)
	t148 = sin(pkin(6));
	t161 = t148 * cos(pkin(7));
	t151 = sin(qJ(2));
	t152 = sin(qJ(1));
	t160 = t151 * t152;
	t154 = cos(qJ(1));
	t159 = t151 * t154;
	t153 = cos(qJ(2));
	t158 = t152 * t153;
	t157 = t153 * t154;
	t156 = qJD(1) * t148;
	t147 = sin(pkin(7));
	t155 = qJD(2) * t147;
	t150 = cos(pkin(6));
	t1 = [0, t154 * t156, -(t150 * t160 - t157) * t155 + (-(-t150 * t157 + t160) * t147 + t154 * t161) * qJD(1), 0, 0, 0; 0, t152 * t156, -(-t150 * t159 - t158) * t155 + (-(-t150 * t158 - t159) * t147 + t152 * t161) * qJD(1), 0, 0, 0; 0, 0, t148 * t151 * t155, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:08
	% EndTime: 2019-10-10 12:08:08
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (55->25), mult. (190->56), div. (0->0), fcn. (202->12), ass. (0->37)
	t235 = cos(pkin(7));
	t231 = sin(pkin(13));
	t234 = cos(pkin(13));
	t237 = sin(qJ(3));
	t240 = cos(qJ(3));
	t230 = -t240 * t231 - t237 * t234;
	t243 = qJD(3) * t230;
	t226 = t235 * t243;
	t259 = qJD(2) * t230 + t226;
	t233 = sin(pkin(6));
	t239 = sin(qJ(1));
	t258 = t233 * t239;
	t242 = cos(qJ(1));
	t257 = t233 * t242;
	t238 = sin(qJ(2));
	t256 = t239 * t238;
	t241 = cos(qJ(2));
	t255 = t239 * t241;
	t254 = t242 * t238;
	t253 = t242 * t241;
	t252 = qJD(1) * t233;
	t250 = t239 * t252;
	t249 = t242 * t252;
	t248 = t231 * t237 - t234 * t240;
	t236 = cos(pkin(6));
	t247 = t236 * t253 - t256;
	t246 = t236 * t255 + t254;
	t245 = t236 * t254 + t255;
	t244 = t236 * t256 - t253;
	t232 = sin(pkin(7));
	t229 = t248 * qJD(3);
	t228 = t248 * t235;
	t227 = t248 * t232;
	t225 = t232 * t243;
	t224 = -t246 * qJD(1) - t245 * qJD(2);
	t223 = -t247 * qJD(1) + t244 * qJD(2);
	t1 = [0, t249, -t223 * t232 + t235 * t249, 0, t244 * t229 + t223 * t228 - t225 * t258 + (t227 * t257 + t245 * t230) * qJD(1) + t259 * t246, 0; 0, t250, -t224 * t232 + t235 * t250, 0, -t245 * t229 + t224 * t228 + t225 * t257 + (t227 * t258 + t244 * t230) * qJD(1) - t259 * t247, 0; 0, 0, t233 * qJD(2) * t238 * t232, 0, -t236 * t225 + (-t226 * t241 - t229 * t238 + (-t228 * t238 - t230 * t241) * qJD(2)) * t233, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:09
	% EndTime: 2019-10-10 12:08:10
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (143->58), mult. (464->118), div. (0->0), fcn. (512->14), ass. (0->48)
	t327 = sin(pkin(13));
	t330 = cos(pkin(13));
	t334 = sin(qJ(3));
	t338 = cos(qJ(3));
	t343 = t327 * t334 - t330 * t338;
	t356 = t343 * qJD(3);
	t328 = sin(pkin(7));
	t335 = sin(qJ(2));
	t355 = t328 * t335;
	t329 = sin(pkin(6));
	t336 = sin(qJ(1));
	t354 = t329 * t336;
	t340 = cos(qJ(1));
	t353 = t329 * t340;
	t352 = t335 * t340;
	t351 = t336 * t335;
	t339 = cos(qJ(2));
	t350 = t336 * t339;
	t349 = t340 * t339;
	t348 = qJD(1) * t336;
	t347 = qJD(1) * t340;
	t346 = t329 * t348;
	t345 = t329 * t347;
	t344 = t327 * t338 + t330 * t334;
	t332 = cos(pkin(6));
	t321 = t332 * t349 - t351;
	t323 = -t332 * t350 - t352;
	t322 = t332 * t352 + t350;
	t342 = t332 * t351 - t349;
	t320 = t344 * qJD(3);
	t337 = cos(qJ(5));
	t333 = sin(qJ(5));
	t331 = cos(pkin(7));
	t318 = t344 * t331;
	t317 = t343 * t331;
	t316 = t344 * t328;
	t315 = t343 * t328;
	t314 = t331 * t356;
	t313 = t331 * t320;
	t312 = t328 * t356;
	t311 = t328 * t320;
	t310 = -t342 * qJD(1) + t321 * qJD(2);
	t309 = t323 * qJD(1) - t322 * qJD(2);
	t308 = -t322 * qJD(1) + t323 * qJD(2);
	t307 = -t321 * qJD(1) + t342 * qJD(2);
	t306 = -t309 * t328 + t331 * t346;
	t305 = -t307 * t328 + t331 * t345;
	t1 = [0, t345, t305, 0, t307 * t317 + t308 * t344 + t323 * t313 + t342 * t356 + (t311 * t336 + t315 * t347) * t329, (t307 * t318 - t308 * t343 - t323 * t314 + t342 * t320 + (-t312 * t336 + t316 * t347) * t329) * t333 - t305 * t337 + ((t316 * t354 + t318 * t323 + t342 * t343) * t337 + (-t323 * t328 + t331 * t354) * t333) * qJD(5); 0, t346, t306, 0, t309 * t317 + t310 * t344 + t321 * t313 - t322 * t356 + (-t311 * t340 + t315 * t348) * t329, (t309 * t318 - t310 * t343 - t321 * t314 - t322 * t320 + (t312 * t340 + t316 * t348) * t329) * t333 - t306 * t337 + ((-t316 * t353 + t321 * t318 - t322 * t343) * t337 + (-t321 * t328 - t331 * t353) * t333) * qJD(5); 0, 0, t329 * qJD(2) * t355, 0, t311 * t332 + (t313 * t339 - t356 * t335 + (-t317 * t335 + t339 * t344) * qJD(2)) * t329, (-t312 * t333 + (t316 * t337 + t331 * t333) * qJD(5)) * t332 + ((-t314 * t339 - t320 * t335) * t333 + ((t318 * t339 - t335 * t343) * t337 - t328 * t339 * t333) * qJD(5) + ((-t318 * t335 - t339 * t343) * t333 - t337 * t355) * qJD(2)) * t329;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end