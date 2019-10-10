% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR12_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR12_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:59
	% EndTime: 2019-10-10 12:49:59
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
	% StartTime: 2019-10-10 12:49:59
	% EndTime: 2019-10-10 12:49:59
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
	% StartTime: 2019-10-10 12:49:59
	% EndTime: 2019-10-10 12:49:59
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (38->24), mult. (136->54), div. (0->0), fcn. (140->10), ass. (0->31)
	t224 = sin(pkin(7));
	t225 = sin(pkin(6));
	t251 = t225 * t224;
	t226 = cos(pkin(7));
	t231 = cos(qJ(3));
	t250 = t226 * t231;
	t228 = sin(qJ(3));
	t232 = cos(qJ(2));
	t249 = t228 * t232;
	t229 = sin(qJ(2));
	t248 = t229 * t231;
	t230 = sin(qJ(1));
	t247 = t230 * t229;
	t246 = t230 * t232;
	t233 = cos(qJ(1));
	t245 = t233 * t229;
	t244 = t233 * t232;
	t243 = qJD(1) * t225;
	t242 = qJD(2) * t228;
	t241 = t230 * t251;
	t240 = t233 * t251;
	t239 = t230 * t243;
	t238 = t233 * t243;
	t227 = cos(pkin(6));
	t237 = t227 * t244 - t247;
	t236 = -t227 * t246 - t245;
	t235 = t227 * t245 + t246;
	t234 = t227 * t247 - t244;
	t223 = t236 * qJD(1) - t235 * qJD(2);
	t222 = -t237 * qJD(1) + t234 * qJD(2);
	t1 = [0, t238, -t222 * t224 + t226 * t238, -t222 * t250 + t236 * t242 + (-t235 * t228 - t231 * t240) * qJD(1) + (-t234 * t231 + (t236 * t226 + t241) * t228) * qJD(3), 0, 0; 0, t239, -t223 * t224 + t226 * t239, -t223 * t250 + t237 * t242 + (-t234 * t228 - t231 * t241) * qJD(1) + (t235 * t231 + (t237 * t226 - t240) * t228) * qJD(3), 0, 0; 0, 0, qJD(2) * t229 * t251, t227 * t224 * qJD(3) * t228 + ((t226 * t249 + t248) * qJD(3) + (t226 * t248 + t249) * qJD(2)) * t225, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:50:00
	% EndTime: 2019-10-10 12:50:01
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (112->50), mult. (332->104), div. (0->0), fcn. (355->12), ass. (0->48)
	t315 = sin(pkin(7));
	t320 = sin(qJ(2));
	t349 = t315 * t320;
	t316 = sin(pkin(6));
	t321 = sin(qJ(1));
	t348 = t316 * t321;
	t324 = cos(qJ(1));
	t347 = t316 * t324;
	t317 = cos(pkin(7));
	t319 = sin(qJ(3));
	t346 = t317 * t319;
	t345 = t319 * t320;
	t323 = cos(qJ(2));
	t344 = t319 * t323;
	t322 = cos(qJ(3));
	t343 = t320 * t322;
	t342 = t321 * t320;
	t341 = t321 * t323;
	t340 = t322 * t323;
	t339 = t324 * t320;
	t338 = t324 * t323;
	t337 = qJD(1) * t316;
	t314 = qJ(4) + pkin(13);
	t312 = sin(t314);
	t336 = qJD(3) * t312;
	t335 = qJD(3) * t315;
	t334 = t321 * t337;
	t333 = t324 * t337;
	t332 = t315 * t334;
	t331 = t315 * t333;
	t318 = cos(pkin(6));
	t308 = t318 * t338 - t342;
	t330 = t308 * t317 - t315 * t347;
	t310 = -t318 * t341 - t339;
	t329 = t310 * t317 + t315 * t348;
	t328 = t317 * t344 + t343;
	t309 = t318 * t339 + t341;
	t327 = t318 * t342 - t338;
	t326 = t309 * t322 + t330 * t319;
	t325 = t329 * t319 - t322 * t327;
	t313 = cos(t314);
	t307 = -t327 * qJD(1) + t308 * qJD(2);
	t306 = t310 * qJD(1) - t309 * qJD(2);
	t305 = -t309 * qJD(1) + t310 * qJD(2);
	t304 = -t308 * qJD(1) + t327 * qJD(2);
	t303 = -t306 * t315 + t317 * t334;
	t302 = -t304 * t315 + t317 * t333;
	t1 = [0, t333, t302, t305 * t319 + (-t304 * t317 - t331) * t322 + t325 * qJD(3), 0, (t304 * t346 + t305 * t322 + t319 * t331) * t312 - t302 * t313 + (t325 * t313 + (-t310 * t315 + t317 * t348) * t312) * qJD(4) + (t319 * t327 + t329 * t322) * t336; 0, t334, t303, t307 * t319 + (-t306 * t317 - t332) * t322 + t326 * qJD(3), 0, (t306 * t346 + t307 * t322 + t319 * t332) * t312 - t303 * t313 + (t326 * t313 + (-t308 * t315 - t317 * t347) * t312) * qJD(4) + (-t309 * t319 + t330 * t322) * t336; 0, 0, t316 * qJD(2) * t349, t318 * t319 * t335 + (t328 * qJD(3) + (t317 * t343 + t344) * qJD(2)) * t316, 0, (t322 * t312 * t335 + (t313 * t315 * t319 + t312 * t317) * qJD(4)) * t318 + ((-t315 * t323 * t312 + t328 * t313) * qJD(4) + (t317 * t340 - t345) * t336 + ((-t317 * t345 + t340) * t312 - t313 * t349) * qJD(2)) * t316;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end