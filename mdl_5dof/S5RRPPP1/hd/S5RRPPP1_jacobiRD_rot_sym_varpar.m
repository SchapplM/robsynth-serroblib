% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPPP1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPPP1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:41
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:41
	% EndTime: 2019-12-29 18:08:42
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:41
	% EndTime: 2019-12-29 18:08:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t34 = sin(qJ(1));
	t41 = qJD(1) * t34;
	t36 = cos(qJ(1));
	t40 = qJD(1) * t36;
	t33 = sin(qJ(2));
	t39 = qJD(2) * t33;
	t35 = cos(qJ(2));
	t38 = qJD(2) * t35;
	t37 = qJD(2) * t36;
	t32 = t34 * t39 - t35 * t40;
	t31 = t33 * t40 + t34 * t38;
	t30 = t33 * t37 + t35 * t41;
	t29 = t33 * t41 - t35 * t37;
	t1 = [t32, t29, 0, 0, 0; -t30, -t31, 0, 0, 0; 0, -t39, 0, 0, 0; t31, t30, 0, 0, 0; t29, t32, 0, 0, 0; 0, -t38, 0, 0, 0; -t41, 0, 0, 0, 0; t40, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:44
	% EndTime: 2019-12-29 18:08:44
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (36->26), mult. (164->65), div. (0->0), fcn. (164->8), ass. (0->28)
	t244 = cos(pkin(5));
	t245 = sin(qJ(2));
	t267 = t244 * t245;
	t247 = cos(qJ(2));
	t266 = t244 * t247;
	t246 = sin(qJ(1));
	t265 = t245 * t246;
	t248 = cos(qJ(1));
	t264 = t245 * t248;
	t263 = t246 * t247;
	t262 = t247 * t248;
	t261 = qJD(1) * t248;
	t260 = qJD(2) * t245;
	t242 = sin(pkin(5));
	t259 = qJD(2) * t247 * t242;
	t241 = sin(pkin(8));
	t243 = cos(pkin(8));
	t258 = -t241 * t245 + t243 * t266;
	t257 = t241 * t266 + t243 * t245;
	t256 = -t242 * t246 + t244 * t264;
	t255 = t242 * t248 + t244 * t265;
	t254 = t258 * t246;
	t253 = t257 * t246;
	t252 = qJD(2) * (t241 * t247 + t243 * t267);
	t251 = qJD(2) * (t241 * t267 - t243 * t247);
	t250 = t258 * qJD(2);
	t249 = t257 * qJD(2);
	t1 = [qJD(2) * t253 + (t256 * t241 - t243 * t262) * qJD(1), qJD(1) * t253 + t248 * t251, 0, 0, 0; -t248 * t249 + (t255 * t241 - t243 * t263) * qJD(1), t246 * t251 - t257 * t261, 0, 0, 0; 0, -t249, 0, 0, 0; qJD(2) * t254 + (t241 * t262 + t256 * t243) * qJD(1), qJD(1) * t254 + t248 * t252, 0, 0, 0; -t248 * t250 + (t241 * t263 + t255 * t243) * qJD(1), t246 * t252 - t258 * t261, 0, 0, 0; 0, -t250, 0, 0, 0; -t246 * t259 + (-t242 * t264 - t244 * t246) * qJD(1), (-qJD(1) * t263 - t248 * t260) * t242, 0, 0, 0; t248 * t259 + (-t242 * t265 + t244 * t248) * qJD(1), (-t246 * t260 + t247 * t261) * t242, 0, 0, 0; 0, t259, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:44
	% EndTime: 2019-12-29 18:08:45
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (36->24), mult. (164->65), div. (0->0), fcn. (164->8), ass. (0->28)
	t293 = cos(pkin(5));
	t294 = sin(qJ(2));
	t316 = t293 * t294;
	t296 = cos(qJ(2));
	t315 = t293 * t296;
	t295 = sin(qJ(1));
	t314 = t294 * t295;
	t297 = cos(qJ(1));
	t313 = t294 * t297;
	t312 = t295 * t296;
	t311 = t296 * t297;
	t310 = qJD(1) * t297;
	t309 = qJD(2) * t294;
	t291 = sin(pkin(5));
	t308 = qJD(2) * t296 * t291;
	t290 = sin(pkin(8));
	t292 = cos(pkin(8));
	t307 = -t290 * t294 + t292 * t315;
	t306 = t290 * t315 + t292 * t294;
	t305 = t291 * t295 - t293 * t313;
	t304 = -t291 * t297 - t293 * t314;
	t303 = t307 * t295;
	t302 = t306 * t295;
	t301 = qJD(2) * (-t290 * t296 - t292 * t316);
	t300 = qJD(2) * (-t290 * t316 + t292 * t296);
	t299 = t307 * qJD(2);
	t298 = t306 * qJD(2);
	t1 = [-t295 * t308 + (-t291 * t313 - t293 * t295) * qJD(1), (-qJD(1) * t312 - t297 * t309) * t291, 0, 0, 0; t297 * t308 + (-t291 * t314 + t293 * t297) * qJD(1), (-t295 * t309 + t296 * t310) * t291, 0, 0, 0; 0, t308, 0, 0, 0; -qJD(2) * t302 + (t305 * t290 + t292 * t311) * qJD(1), -qJD(1) * t302 + t297 * t300, 0, 0, 0; t297 * t298 + (t304 * t290 + t292 * t312) * qJD(1), t295 * t300 + t306 * t310, 0, 0, 0; 0, t298, 0, 0, 0; -qJD(2) * t303 + (-t290 * t311 + t305 * t292) * qJD(1), -qJD(1) * t303 + t297 * t301, 0, 0, 0; t297 * t299 + (-t290 * t312 + t304 * t292) * qJD(1), t295 * t301 + t307 * t310, 0, 0, 0; 0, t299, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:08:45
	% EndTime: 2019-12-29 18:08:45
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (36->25), mult. (164->65), div. (0->0), fcn. (164->8), ass. (0->28)
	t330 = cos(pkin(5));
	t331 = sin(qJ(2));
	t353 = t330 * t331;
	t333 = cos(qJ(2));
	t352 = t330 * t333;
	t332 = sin(qJ(1));
	t351 = t331 * t332;
	t334 = cos(qJ(1));
	t350 = t331 * t334;
	t349 = t332 * t333;
	t348 = t333 * t334;
	t347 = qJD(1) * t334;
	t346 = qJD(2) * t331;
	t328 = sin(pkin(5));
	t345 = qJD(2) * t333 * t328;
	t327 = sin(pkin(8));
	t329 = cos(pkin(8));
	t344 = -t327 * t331 + t329 * t352;
	t343 = t327 * t352 + t329 * t331;
	t342 = -t328 * t332 + t330 * t350;
	t341 = t328 * t334 + t330 * t351;
	t340 = t344 * t332;
	t339 = t343 * t332;
	t338 = qJD(2) * (-t327 * t333 - t329 * t353);
	t337 = qJD(2) * (t327 * t353 - t329 * t333);
	t336 = t344 * qJD(2);
	t335 = t343 * qJD(2);
	t1 = [-t332 * t345 + (-t328 * t350 - t330 * t332) * qJD(1), (-qJD(1) * t349 - t334 * t346) * t328, 0, 0, 0; t334 * t345 + (-t328 * t351 + t330 * t334) * qJD(1), (-t332 * t346 + t333 * t347) * t328, 0, 0, 0; 0, t345, 0, 0, 0; -qJD(2) * t340 + (-t327 * t348 - t342 * t329) * qJD(1), -qJD(1) * t340 + t334 * t338, 0, 0, 0; t334 * t336 + (-t327 * t349 - t341 * t329) * qJD(1), t332 * t338 + t344 * t347, 0, 0, 0; 0, t336, 0, 0, 0; qJD(2) * t339 + (t342 * t327 - t329 * t348) * qJD(1), qJD(1) * t339 + t334 * t337, 0, 0, 0; -t334 * t335 + (t341 * t327 - t329 * t349) * qJD(1), t332 * t337 - t343 * t347, 0, 0, 0; 0, -t335, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end