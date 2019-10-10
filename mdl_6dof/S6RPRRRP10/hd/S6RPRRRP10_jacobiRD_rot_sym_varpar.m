% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRRP10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:12
	% EndTime: 2019-10-10 08:54:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:12
	% EndTime: 2019-10-10 08:54:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0, 0; -t31, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t31, 0, 0, 0, 0, 0; -t30, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:12
	% EndTime: 2019-10-10 08:54:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:13
	% EndTime: 2019-10-10 08:54:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (11->9), mult. (36->13), div. (0->0), fcn. (36->4), ass. (0->14)
	t35 = sin(qJ(1));
	t42 = qJD(1) * t35;
	t37 = cos(qJ(1));
	t41 = qJD(1) * t37;
	t34 = sin(qJ(3));
	t40 = qJD(3) * t34;
	t36 = cos(qJ(3));
	t39 = qJD(3) * t36;
	t38 = qJD(3) * t37;
	t33 = -t35 * t40 + t36 * t41;
	t32 = t34 * t41 + t35 * t39;
	t31 = t34 * t38 + t36 * t42;
	t30 = -t34 * t42 + t36 * t38;
	t1 = [t30, 0, t33, 0, 0, 0; t32, 0, t31, 0, 0, 0; 0, 0, -t39, 0, 0, 0; -t31, 0, -t32, 0, 0, 0; t33, 0, t30, 0, 0, 0; 0, 0, t40, 0, 0, 0; -t41, 0, 0, 0, 0, 0; -t42, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:14
	% EndTime: 2019-10-10 08:54:14
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (49->24), mult. (173->45), div. (0->0), fcn. (173->6), ass. (0->32)
	t239 = cos(qJ(1));
	t235 = sin(qJ(3));
	t244 = qJD(1) * t235 + qJD(4);
	t260 = t244 * t239;
	t236 = sin(qJ(1));
	t238 = cos(qJ(3));
	t254 = qJD(3) * t239;
	t246 = t238 * t254;
	t259 = t244 * t236 - t246;
	t258 = qJD(1) * t236;
	t257 = qJD(1) * t239;
	t256 = qJD(3) * t235;
	t255 = qJD(3) * t238;
	t253 = qJD(4) * t235;
	t252 = qJD(4) * t238;
	t251 = t238 * t257;
	t237 = cos(qJ(4));
	t250 = t237 * t255;
	t249 = t236 * t255;
	t234 = sin(qJ(4));
	t248 = t234 * t252;
	t247 = t237 * t252;
	t245 = -qJD(1) - t253;
	t243 = t245 * t239;
	t242 = t236 * t256 - t251;
	t241 = t235 * t254 + t238 * t258;
	t240 = t237 * t256 + t248;
	t233 = t237 * t260 + (t245 * t234 + t250) * t236;
	t232 = t245 * t237 * t236 + (-t249 - t260) * t234;
	t231 = t234 * t243 - t259 * t237;
	t230 = t259 * t234 + t237 * t243;
	t1 = [t231, 0, -t240 * t236 + t237 * t251, t232, 0, 0; t233, 0, t241 * t237 + t239 * t248, -t230, 0, 0; 0, 0, t234 * t253 - t250, t234 * t256 - t247, 0, 0; t230, 0, t242 * t234 - t236 * t247, -t233, 0, 0; t232, 0, -t241 * t234 + t239 * t247, t231, 0, 0; 0, 0, t234 * t255 + t237 * t253, t240, 0, 0; t241, 0, t235 * t257 + t249, 0, 0, 0; t242, 0, t235 * t258 - t246, 0, 0, 0; 0, 0, -t256, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:14
	% EndTime: 2019-10-10 08:54:14
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (185->30), mult. (233->47), div. (0->0), fcn. (233->6), ass. (0->34)
	t297 = sin(qJ(1));
	t294 = qJD(4) + qJD(5);
	t296 = sin(qJ(3));
	t302 = qJD(1) * t296 + t294;
	t298 = cos(qJ(3));
	t299 = cos(qJ(1));
	t310 = qJD(3) * t299;
	t304 = t298 * t310;
	t318 = t302 * t297 - t304;
	t311 = qJD(3) * t298;
	t305 = t297 * t311;
	t317 = t302 * t299 + t305;
	t316 = t294 * t296;
	t315 = t294 * t298;
	t314 = qJD(1) * t297;
	t313 = qJD(1) * t299;
	t312 = qJD(3) * t296;
	t295 = qJ(4) + qJ(5);
	t292 = sin(t295);
	t309 = t292 * t316;
	t308 = t292 * t315;
	t293 = cos(t295);
	t307 = t293 * t315;
	t306 = t299 * t294 * t293;
	t303 = -qJD(1) - t316;
	t301 = -t297 * t312 + t298 * t313;
	t300 = t296 * t310 + t298 * t314;
	t286 = t293 * t312 + t308;
	t285 = t292 * t312 - t307;
	t284 = -t292 * t314 + t317 * t293 - t297 * t309;
	t283 = t303 * t297 * t293 - t317 * t292;
	t282 = t303 * t299 * t292 - t318 * t293;
	t281 = t318 * t292 - t293 * t313 - t296 * t306;
	t1 = [t282, 0, t301 * t293 - t297 * t308, t283, t283, 0; t284, 0, t300 * t293 + t299 * t308, -t281, -t281, 0; 0, 0, -t293 * t311 + t309, t285, t285, 0; t281, 0, -t301 * t292 - t297 * t307, -t284, -t284, 0; t283, 0, -t300 * t292 + t298 * t306, t282, t282, 0; 0, 0, t292 * t311 + t293 * t316, t286, t286, 0; t300, 0, t296 * t313 + t305, 0, 0, 0; -t301, 0, t296 * t314 - t304, 0, 0, 0; 0, 0, -t312, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:54:15
	% EndTime: 2019-10-10 08:54:15
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (185->27), mult. (233->44), div. (0->0), fcn. (233->6), ass. (0->34)
	t353 = sin(qJ(1));
	t350 = qJD(4) + qJD(5);
	t352 = sin(qJ(3));
	t372 = t350 * t352;
	t361 = qJD(1) + t372;
	t374 = t353 * t361;
	t360 = qJD(1) * t352 + t350;
	t354 = cos(qJ(3));
	t355 = cos(qJ(1));
	t366 = qJD(3) * t355;
	t362 = t354 * t366;
	t373 = t360 * t353 - t362;
	t371 = t350 * t354;
	t370 = qJD(1) * t353;
	t369 = qJD(1) * t355;
	t368 = qJD(3) * t352;
	t367 = qJD(3) * t354;
	t351 = qJ(4) + qJ(5);
	t348 = sin(t351);
	t365 = t348 * t371;
	t349 = cos(t351);
	t364 = t349 * t371;
	t363 = t353 * t367;
	t359 = t361 * t355;
	t358 = -t353 * t368 + t354 * t369;
	t357 = t352 * t366 + t354 * t370;
	t356 = t360 * t355 + t363;
	t347 = -t349 * t368 - t365;
	t346 = t348 * t368 - t364;
	t345 = -t348 * t374 + t356 * t349;
	t344 = t356 * t348 + t349 * t374;
	t343 = t348 * t359 + t373 * t349;
	t342 = -t373 * t348 + t349 * t359;
	t1 = [-t343, 0, t358 * t349 - t353 * t365, -t344, -t344, 0; t345, 0, t357 * t349 + t355 * t365, t342, t342, 0; 0, 0, t348 * t372 - t349 * t367, t346, t346, 0; t357, 0, t352 * t369 + t363, 0, 0, 0; -t358, 0, t352 * t370 - t362, 0, 0, 0; 0, 0, -t368, 0, 0, 0; t342, 0, t358 * t348 + t353 * t364, t345, t345, 0; t344, 0, t357 * t348 - t355 * t364, t343, t343, 0; 0, 0, -t348 * t367 - t349 * t372, t347, t347, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end