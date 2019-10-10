% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR13
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
% Datum: 2019-10-10 12:14
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRR13_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR13_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:18
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
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:18
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (8->8), mult. (35->23), div. (0->0), fcn. (35->8), ass. (0->15)
	t192 = sin(pkin(6));
	t205 = t192 * cos(pkin(7));
	t195 = sin(qJ(2));
	t196 = sin(qJ(1));
	t204 = t195 * t196;
	t198 = cos(qJ(1));
	t203 = t195 * t198;
	t197 = cos(qJ(2));
	t202 = t196 * t197;
	t201 = t197 * t198;
	t200 = qJD(1) * t192;
	t191 = sin(pkin(7));
	t199 = qJD(2) * t191;
	t194 = cos(pkin(6));
	t1 = [0, t198 * t200, -(t194 * t204 - t201) * t199 + (-(-t194 * t201 + t204) * t191 + t198 * t205) * qJD(1), 0, 0, 0; 0, t196 * t200, -(-t194 * t203 - t202) * t199 + (-(-t194 * t202 - t203) * t191 + t196 * t205) * qJD(1), 0, 0, 0; 0, 0, t192 * t195 * t199, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:18
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (38->24), mult. (136->54), div. (0->0), fcn. (140->10), ass. (0->31)
	t223 = sin(pkin(7));
	t224 = sin(pkin(6));
	t250 = t224 * t223;
	t225 = cos(pkin(7));
	t230 = cos(qJ(3));
	t249 = t225 * t230;
	t227 = sin(qJ(3));
	t231 = cos(qJ(2));
	t248 = t227 * t231;
	t228 = sin(qJ(2));
	t247 = t228 * t230;
	t229 = sin(qJ(1));
	t246 = t229 * t228;
	t245 = t229 * t231;
	t232 = cos(qJ(1));
	t244 = t232 * t228;
	t243 = t232 * t231;
	t242 = qJD(1) * t224;
	t241 = qJD(2) * t227;
	t240 = t229 * t250;
	t239 = t232 * t250;
	t238 = t229 * t242;
	t237 = t232 * t242;
	t226 = cos(pkin(6));
	t236 = t226 * t243 - t246;
	t235 = -t226 * t245 - t244;
	t234 = t226 * t244 + t245;
	t233 = t226 * t246 - t243;
	t222 = t235 * qJD(1) - t234 * qJD(2);
	t221 = -t236 * qJD(1) + t233 * qJD(2);
	t1 = [0, t237, -t221 * t223 + t225 * t237, 0, -t221 * t249 + t235 * t241 + (-t234 * t227 - t230 * t239) * qJD(1) + (-t233 * t230 + (t235 * t225 + t240) * t227) * qJD(3), 0; 0, t238, -t222 * t223 + t225 * t238, 0, -t222 * t249 + t236 * t241 + (-t233 * t227 - t230 * t240) * qJD(1) + (t234 * t230 + (t236 * t225 - t239) * t227) * qJD(3), 0; 0, 0, qJD(2) * t228 * t250, 0, t226 * t223 * qJD(3) * t227 + ((t225 * t248 + t247) * qJD(3) + (t225 * t247 + t248) * qJD(2)) * t224, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:20
	% EndTime: 2019-10-10 12:14:20
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (112->50), mult. (332->104), div. (0->0), fcn. (355->12), ass. (0->48)
	t314 = sin(pkin(7));
	t319 = sin(qJ(2));
	t348 = t314 * t319;
	t315 = sin(pkin(6));
	t320 = sin(qJ(1));
	t347 = t315 * t320;
	t323 = cos(qJ(1));
	t346 = t315 * t323;
	t316 = cos(pkin(7));
	t318 = sin(qJ(3));
	t345 = t316 * t318;
	t344 = t318 * t319;
	t322 = cos(qJ(2));
	t343 = t318 * t322;
	t321 = cos(qJ(3));
	t342 = t319 * t321;
	t341 = t320 * t319;
	t340 = t320 * t322;
	t339 = t321 * t322;
	t338 = t323 * t319;
	t337 = t323 * t322;
	t336 = qJD(1) * t315;
	t313 = pkin(13) + qJ(5);
	t311 = sin(t313);
	t335 = qJD(3) * t311;
	t334 = qJD(3) * t314;
	t333 = t320 * t336;
	t332 = t323 * t336;
	t331 = t314 * t333;
	t330 = t314 * t332;
	t317 = cos(pkin(6));
	t307 = t317 * t337 - t341;
	t329 = t307 * t316 - t314 * t346;
	t309 = -t317 * t340 - t338;
	t328 = t309 * t316 + t314 * t347;
	t327 = t316 * t343 + t342;
	t308 = t317 * t338 + t340;
	t326 = t317 * t341 - t337;
	t325 = t308 * t321 + t329 * t318;
	t324 = t328 * t318 - t321 * t326;
	t312 = cos(t313);
	t306 = -t326 * qJD(1) + t307 * qJD(2);
	t305 = t309 * qJD(1) - t308 * qJD(2);
	t304 = -t308 * qJD(1) + t309 * qJD(2);
	t303 = -t307 * qJD(1) + t326 * qJD(2);
	t302 = -t305 * t314 + t316 * t333;
	t301 = -t303 * t314 + t316 * t332;
	t1 = [0, t332, t301, 0, t304 * t318 + (-t303 * t316 - t330) * t321 + t324 * qJD(3), (t303 * t345 + t304 * t321 + t318 * t330) * t311 - t301 * t312 + (t324 * t312 + (-t309 * t314 + t316 * t347) * t311) * qJD(5) + (t318 * t326 + t328 * t321) * t335; 0, t333, t302, 0, t306 * t318 + (-t305 * t316 - t331) * t321 + t325 * qJD(3), (t305 * t345 + t306 * t321 + t318 * t331) * t311 - t302 * t312 + (t325 * t312 + (-t307 * t314 - t316 * t346) * t311) * qJD(5) + (-t308 * t318 + t329 * t321) * t335; 0, 0, t315 * qJD(2) * t348, 0, t317 * t318 * t334 + (t327 * qJD(3) + (t316 * t342 + t343) * qJD(2)) * t315, (t321 * t311 * t334 + (t312 * t314 * t318 + t311 * t316) * qJD(5)) * t317 + ((-t314 * t322 * t311 + t327 * t312) * qJD(5) + (t316 * t339 - t344) * t335 + ((-t316 * t344 + t339) * t311 - t312 * t348) * qJD(2)) * t315;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end