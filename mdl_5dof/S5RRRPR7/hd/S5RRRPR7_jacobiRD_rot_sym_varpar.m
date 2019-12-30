% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRPR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:18
	% EndTime: 2019-12-29 20:05:18
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.06s
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
	% StartTime: 2019-12-29 20:05:13
	% EndTime: 2019-12-29 20:05:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (61->14), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t71 = qJ(2) + qJ(3);
	t68 = sin(t71);
	t70 = qJD(2) + qJD(3);
	t79 = t70 * t68;
	t69 = cos(t71);
	t78 = t70 * t69;
	t72 = sin(qJ(1));
	t77 = t70 * t72;
	t73 = cos(qJ(1));
	t76 = t70 * t73;
	t75 = qJD(1) * t72;
	t74 = qJD(1) * t73;
	t67 = t68 * t77 - t69 * t74;
	t66 = t68 * t74 + t69 * t77;
	t65 = t68 * t76 + t69 * t75;
	t64 = t68 * t75 - t69 * t76;
	t1 = [t67, t64, t64, 0, 0; -t65, -t66, -t66, 0, 0; 0, -t79, -t79, 0, 0; t66, t65, t65, 0, 0; t64, t67, t67, 0, 0; 0, -t78, -t78, 0, 0; -t75, 0, 0, 0, 0; t74, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:15
	% EndTime: 2019-12-29 20:05:15
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (90->22), mult. (114->36), div. (0->0), fcn. (114->6), ass. (0->31)
	t282 = qJ(2) + qJ(3);
	t279 = sin(t282);
	t281 = qJD(2) + qJD(3);
	t300 = t279 * t281;
	t285 = sin(qJ(1));
	t299 = t281 * t285;
	t286 = cos(qJ(1));
	t298 = t281 * t286;
	t283 = sin(pkin(9));
	t297 = t283 * t285;
	t296 = t283 * t286;
	t284 = cos(pkin(9));
	t295 = t284 * t285;
	t294 = t284 * t286;
	t293 = qJD(1) * t285;
	t292 = qJD(1) * t286;
	t291 = t284 * t300;
	t290 = t279 * t299;
	t289 = t279 * t298;
	t280 = cos(t282);
	t288 = t279 * t292 + t280 * t299;
	t287 = t279 * t293 - t280 * t298;
	t278 = t281 * t280;
	t277 = t283 * t300;
	t276 = t280 * t292 - t290;
	t275 = -t280 * t293 - t289;
	t274 = t288 * t284;
	t273 = t288 * t283;
	t272 = t287 * t284;
	t271 = t287 * t283;
	t1 = [t284 * t290 + (-t280 * t294 - t297) * qJD(1), t272, t272, 0, 0; -t284 * t289 + (-t280 * t295 + t296) * qJD(1), -t274, -t274, 0, 0; 0, -t291, -t291, 0, 0; -t283 * t290 + (t280 * t296 - t295) * qJD(1), -t271, -t271, 0, 0; t283 * t289 + (t280 * t297 + t294) * qJD(1), t273, t273, 0, 0; 0, t277, t277, 0, 0; -t288, t275, t275, 0, 0; -t287, t276, t276, 0, 0; 0, t278, t278, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:05:16
	% EndTime: 2019-12-29 20:05:16
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (240->27), mult. (226->45), div. (0->0), fcn. (226->6), ass. (0->42)
	t343 = cos(qJ(1));
	t341 = qJ(2) + qJ(3);
	t338 = cos(t341);
	t354 = qJD(5) * t338;
	t349 = -qJD(1) + t354;
	t362 = t343 * t349;
	t348 = qJD(1) * t338 - qJD(5);
	t337 = sin(t341);
	t340 = qJD(2) + qJD(3);
	t342 = sin(qJ(1));
	t359 = t340 * t342;
	t353 = t337 * t359;
	t361 = t348 * t343 - t353;
	t360 = t337 * t340;
	t334 = t340 * t338;
	t358 = t340 * t343;
	t357 = qJD(1) * t342;
	t356 = qJD(1) * t343;
	t355 = qJD(5) * t337;
	t352 = t337 * t358;
	t339 = pkin(9) + qJ(5);
	t335 = sin(t339);
	t351 = t335 * t355;
	t336 = cos(t339);
	t350 = t336 * t355;
	t347 = t349 * t342;
	t346 = t337 * t356 + t338 * t359;
	t345 = t337 * t357 - t338 * t358;
	t344 = t348 * t342 + t352;
	t333 = t338 * t356 - t353;
	t332 = -t338 * t357 - t352;
	t331 = -t335 * t354 - t336 * t360;
	t330 = t335 * t360 - t336 * t354;
	t329 = -t346 * t336 + t342 * t351;
	t328 = t346 * t335 + t342 * t350;
	t327 = t345 * t336 + t343 * t351;
	t326 = -t345 * t335 + t343 * t350;
	t325 = t335 * t347 - t361 * t336;
	t324 = t361 * t335 + t336 * t347;
	t323 = t335 * t362 + t344 * t336;
	t322 = t344 * t335 - t336 * t362;
	t1 = [t325, t327, t327, 0, t322; -t323, t329, t329, 0, -t324; 0, t331, t331, 0, -t335 * t334 - t350; t324, t326, t326, 0, t323; t322, t328, t328, 0, t325; 0, t330, t330, 0, -t336 * t334 + t351; -t346, t332, t332, 0, 0; -t345, t333, t333, 0, 0; 0, t334, t334, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end