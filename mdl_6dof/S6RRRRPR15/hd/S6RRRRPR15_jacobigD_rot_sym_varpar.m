% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:56
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR15_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR15_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
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
	% StartTime: 2019-10-10 12:56:27
	% EndTime: 2019-10-10 12:56:27
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (38->24), mult. (136->54), div. (0->0), fcn. (140->10), ass. (0->31)
	t215 = sin(pkin(7));
	t216 = sin(pkin(6));
	t242 = t216 * t215;
	t217 = cos(pkin(7));
	t222 = cos(qJ(3));
	t241 = t217 * t222;
	t219 = sin(qJ(3));
	t223 = cos(qJ(2));
	t240 = t219 * t223;
	t220 = sin(qJ(2));
	t239 = t220 * t222;
	t221 = sin(qJ(1));
	t238 = t221 * t220;
	t237 = t221 * t223;
	t224 = cos(qJ(1));
	t236 = t224 * t220;
	t235 = t224 * t223;
	t234 = qJD(1) * t216;
	t233 = qJD(2) * t219;
	t232 = t221 * t242;
	t231 = t224 * t242;
	t230 = t221 * t234;
	t229 = t224 * t234;
	t218 = cos(pkin(6));
	t228 = t218 * t235 - t238;
	t227 = -t218 * t237 - t236;
	t226 = t218 * t236 + t237;
	t225 = t218 * t238 - t235;
	t214 = t227 * qJD(1) - t226 * qJD(2);
	t213 = -t228 * qJD(1) + t225 * qJD(2);
	t1 = [0, t229, -t213 * t215 + t217 * t229, -t213 * t241 + t227 * t233 + (-t226 * t219 - t222 * t231) * qJD(1) + (-t225 * t222 + (t227 * t217 + t232) * t219) * qJD(3), 0, 0; 0, t230, -t214 * t215 + t217 * t230, -t214 * t241 + t228 * t233 + (-t225 * t219 - t222 * t232) * qJD(1) + (t226 * t222 + (t228 * t217 - t231) * t219) * qJD(3), 0, 0; 0, 0, qJD(2) * t220 * t242, t218 * t215 * qJD(3) * t219 + ((t217 * t240 + t239) * qJD(3) + (t217 * t239 + t240) * qJD(2)) * t216, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:27
	% EndTime: 2019-10-10 12:56:27
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (38->24), mult. (136->54), div. (0->0), fcn. (140->10), ass. (0->31)
	t240 = sin(pkin(7));
	t241 = sin(pkin(6));
	t267 = t241 * t240;
	t242 = cos(pkin(7));
	t247 = cos(qJ(3));
	t266 = t242 * t247;
	t244 = sin(qJ(3));
	t248 = cos(qJ(2));
	t265 = t244 * t248;
	t245 = sin(qJ(2));
	t264 = t245 * t247;
	t246 = sin(qJ(1));
	t263 = t246 * t245;
	t262 = t246 * t248;
	t249 = cos(qJ(1));
	t261 = t249 * t245;
	t260 = t249 * t248;
	t259 = qJD(1) * t241;
	t258 = qJD(2) * t244;
	t257 = t246 * t267;
	t256 = t249 * t267;
	t255 = t246 * t259;
	t254 = t249 * t259;
	t243 = cos(pkin(6));
	t253 = t243 * t260 - t263;
	t252 = -t243 * t262 - t261;
	t251 = t243 * t261 + t262;
	t250 = t243 * t263 - t260;
	t239 = t252 * qJD(1) - t251 * qJD(2);
	t238 = -t253 * qJD(1) + t250 * qJD(2);
	t1 = [0, t254, -t238 * t240 + t242 * t254, -t238 * t266 + t252 * t258 + (-t251 * t244 - t247 * t256) * qJD(1) + (-t250 * t247 + (t252 * t242 + t257) * t244) * qJD(3), 0, 0; 0, t255, -t239 * t240 + t242 * t255, -t239 * t266 + t253 * t258 + (-t250 * t244 - t247 * t257) * qJD(1) + (t251 * t247 + (t253 * t242 - t256) * t244) * qJD(3), 0, 0; 0, 0, qJD(2) * t245 * t267, t243 * t240 * qJD(3) * t244 + ((t242 * t265 + t264) * qJD(3) + (t242 * t264 + t265) * qJD(2)) * t241, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:28
	% EndTime: 2019-10-10 12:56:28
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (100->49), mult. (332->103), div. (0->0), fcn. (355->12), ass. (0->47)
	t297 = sin(pkin(7));
	t302 = sin(qJ(3));
	t333 = t297 * t302;
	t303 = sin(qJ(2));
	t332 = t297 * t303;
	t298 = sin(pkin(6));
	t304 = sin(qJ(1));
	t331 = t298 * t304;
	t308 = cos(qJ(1));
	t330 = t298 * t308;
	t299 = cos(pkin(7));
	t329 = t299 * t302;
	t328 = t302 * t303;
	t307 = cos(qJ(2));
	t327 = t302 * t307;
	t306 = cos(qJ(3));
	t326 = t303 * t306;
	t325 = t304 * t303;
	t324 = t304 * t307;
	t323 = t306 * t307;
	t322 = t308 * t303;
	t321 = t308 * t307;
	t320 = qJD(1) * t298;
	t305 = cos(qJ(4));
	t319 = qJD(3) * t305;
	t318 = t304 * t320;
	t317 = t308 * t320;
	t316 = t297 * t318;
	t315 = t297 * t317;
	t300 = cos(pkin(6));
	t293 = t300 * t321 - t325;
	t314 = t293 * t299 - t297 * t330;
	t295 = -t300 * t324 - t322;
	t313 = t295 * t299 + t297 * t331;
	t312 = t299 * t327 + t326;
	t294 = t300 * t322 + t324;
	t311 = t300 * t325 - t321;
	t310 = t294 * t306 + t314 * t302;
	t309 = t313 * t302 - t306 * t311;
	t301 = sin(qJ(4));
	t292 = -t311 * qJD(1) + t293 * qJD(2);
	t291 = t295 * qJD(1) - t294 * qJD(2);
	t290 = -t294 * qJD(1) + t295 * qJD(2);
	t289 = -t293 * qJD(1) + t311 * qJD(2);
	t288 = -t291 * t297 + t299 * t318;
	t287 = -t289 * t297 + t299 * t317;
	t1 = [0, t317, t287, t290 * t302 + (-t289 * t299 - t315) * t306 + t309 * qJD(3), 0, (t289 * t329 + t290 * t306 + t302 * t315) * t305 + t287 * t301 + (-t309 * t301 + (-t295 * t297 + t299 * t331) * t305) * qJD(4) + (t302 * t311 + t313 * t306) * t319; 0, t318, t288, t292 * t302 + (-t291 * t299 - t316) * t306 + t310 * qJD(3), 0, (t291 * t329 + t292 * t306 + t302 * t316) * t305 + t288 * t301 + (-t310 * t301 + (-t293 * t297 - t299 * t330) * t305) * qJD(4) + (-t294 * t302 + t314 * t306) * t319; 0, 0, t298 * qJD(2) * t332, t300 * qJD(3) * t333 + (t312 * qJD(3) + (t299 * t326 + t327) * qJD(2)) * t298, 0, (t297 * t306 * t319 + (t299 * t305 - t301 * t333) * qJD(4)) * t300 + ((-t297 * t307 * t305 - t312 * t301) * qJD(4) + (t299 * t323 - t328) * t319 + ((-t299 * t328 + t323) * t305 + t301 * t332) * qJD(2)) * t298;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end