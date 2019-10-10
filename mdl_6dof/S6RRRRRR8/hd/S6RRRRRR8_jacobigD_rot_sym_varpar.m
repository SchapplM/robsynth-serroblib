% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR8_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR8_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobigD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:31
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
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:31
	% DurationCPUTime: 0.13s
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
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:31
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (68->24), mult. (237->54), div. (0->0), fcn. (245->10), ass. (0->34)
	t239 = sin(pkin(7));
	t240 = sin(pkin(6));
	t266 = t240 * t239;
	t241 = cos(pkin(7));
	t246 = cos(qJ(3));
	t265 = t241 * t246;
	t243 = sin(qJ(3));
	t247 = cos(qJ(2));
	t264 = t243 * t247;
	t244 = sin(qJ(2));
	t263 = t244 * t246;
	t245 = sin(qJ(1));
	t262 = t245 * t244;
	t261 = t245 * t247;
	t248 = cos(qJ(1));
	t260 = t248 * t244;
	t259 = t248 * t247;
	t258 = qJD(1) * t240;
	t257 = qJD(2) * t243;
	t256 = t245 * t266;
	t255 = t248 * t266;
	t254 = t245 * t258;
	t253 = t248 * t258;
	t242 = cos(pkin(6));
	t252 = t242 * t259 - t262;
	t251 = -t242 * t261 - t260;
	t250 = t242 * t260 + t261;
	t249 = t242 * t262 - t259;
	t238 = t251 * qJD(1) - t250 * qJD(2);
	t237 = -t252 * qJD(1) + t249 * qJD(2);
	t236 = t242 * t239 * qJD(3) * t243 + ((t241 * t264 + t263) * qJD(3) + (t241 * t263 + t264) * qJD(2)) * t240;
	t235 = -t238 * t265 + t252 * t257 + (-t249 * t243 - t246 * t256) * qJD(1) + (t250 * t246 + (t252 * t241 - t255) * t243) * qJD(3);
	t234 = -t237 * t265 + t251 * t257 + (-t250 * t243 - t246 * t255) * qJD(1) + (-t249 * t246 + (t251 * t241 + t256) * t243) * qJD(3);
	t1 = [0, t253, -t237 * t239 + t241 * t253, t234, t234, 0; 0, t254, -t238 * t239 + t241 * t254, t235, t235, 0; 0, 0, qJD(2) * t244 * t266, t236, t236, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:33
	% EndTime: 2019-10-10 13:29:33
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (148->49), mult. (433->99), div. (0->0), fcn. (460->12), ass. (0->51)
	t346 = sin(pkin(7));
	t351 = sin(qJ(2));
	t379 = t346 * t351;
	t347 = sin(pkin(6));
	t352 = sin(qJ(1));
	t378 = t347 * t352;
	t355 = cos(qJ(1));
	t377 = t347 * t355;
	t350 = sin(qJ(3));
	t376 = t350 * t351;
	t354 = cos(qJ(2));
	t375 = t350 * t354;
	t353 = cos(qJ(3));
	t374 = t351 * t353;
	t373 = t352 * t351;
	t372 = t352 * t354;
	t371 = t353 * t354;
	t370 = t355 * t351;
	t369 = t355 * t354;
	t368 = qJD(1) * t347;
	t345 = qJ(4) + qJ(5);
	t342 = sin(t345);
	t367 = qJD(3) * t342;
	t366 = qJD(3) * t346;
	t365 = t352 * t368;
	t364 = t355 * t368;
	t349 = cos(pkin(6));
	t338 = t349 * t369 - t373;
	t348 = cos(pkin(7));
	t363 = t338 * t348 - t346 * t377;
	t340 = -t349 * t372 - t370;
	t362 = t340 * t348 + t346 * t378;
	t361 = t348 * t375 + t374;
	t339 = t349 * t370 + t372;
	t360 = t349 * t373 - t369;
	t334 = -qJD(1) * t338 + qJD(2) * t360;
	t359 = t334 * t348 + t346 * t364;
	t336 = qJD(1) * t340 - qJD(2) * t339;
	t358 = t336 * t348 + t346 * t365;
	t357 = t339 * t353 + t350 * t363;
	t356 = t350 * t362 - t353 * t360;
	t344 = qJD(4) + qJD(5);
	t343 = cos(t345);
	t337 = -qJD(1) * t360 + qJD(2) * t338;
	t335 = -qJD(1) * t339 + qJD(2) * t340;
	t333 = -t336 * t346 + t348 * t365;
	t332 = -t334 * t346 + t348 * t364;
	t331 = t349 * t350 * t366 + (t361 * qJD(3) + (t348 * t374 + t375) * qJD(2)) * t347;
	t330 = qJD(3) * t357 + t337 * t350 - t353 * t358;
	t329 = qJD(3) * t356 + t335 * t350 - t353 * t359;
	t1 = [0, t364, t332, t329, t329, (t344 * t356 - t332) * t343 + (t335 * t353 + (-t340 * t346 + t348 * t378) * t344 + t359 * t350) * t342 + (t350 * t360 + t353 * t362) * t367; 0, t365, t333, t330, t330, (t344 * t357 - t333) * t343 + (t337 * t353 + (-t338 * t346 - t348 * t377) * t344 + t358 * t350) * t342 + (-t339 * t350 + t353 * t363) * t367; 0, 0, t347 * qJD(2) * t379, t331, t331, (t346 * t350 * t344 * t343 + (t344 * t348 + t353 * t366) * t342) * t349 + ((-t346 * t354 * t342 + t343 * t361) * t344 + (t348 * t371 - t376) * t367 + ((-t348 * t376 + t371) * t342 - t343 * t379) * qJD(2)) * t347;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end