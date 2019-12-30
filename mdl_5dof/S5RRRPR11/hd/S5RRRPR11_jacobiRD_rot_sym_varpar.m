% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRPR11_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR11_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:25
	% EndTime: 2019-12-29 20:16:26
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:14
	% EndTime: 2019-12-29 20:16:14
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
	% StartTime: 2019-12-29 20:16:19
	% EndTime: 2019-12-29 20:16:19
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
	% StartTime: 2019-12-29 20:16:21
	% EndTime: 2019-12-29 20:16:21
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t232 = cos(qJ(3));
	t234 = cos(qJ(1));
	t256 = t232 * t234;
	t231 = sin(qJ(1));
	t255 = qJD(1) * t231;
	t233 = cos(qJ(2));
	t254 = qJD(1) * t233;
	t253 = qJD(1) * t234;
	t230 = sin(qJ(2));
	t252 = qJD(2) * t230;
	t251 = qJD(2) * t233;
	t250 = qJD(2) * t234;
	t229 = sin(qJ(3));
	t249 = qJD(3) * t229;
	t248 = qJD(3) * t230;
	t247 = qJD(3) * t233;
	t246 = t232 * t252;
	t245 = t232 * t248;
	t244 = t231 * t252;
	t243 = t231 * t251;
	t242 = t230 * t250;
	t241 = t233 * t250;
	t240 = -qJD(1) + t247;
	t239 = -qJD(3) + t254;
	t238 = t240 * t229;
	t237 = t230 * t253 + t243;
	t236 = -t230 * t255 + t241;
	t235 = t239 * t231 + t242;
	t228 = -t239 * t256 + (t238 + t246) * t231;
	t227 = t240 * t232 * t231 + (t239 * t234 - t244) * t229;
	t226 = t235 * t232 + t234 * t238;
	t225 = t235 * t229 - t240 * t256;
	t1 = [t228, -t232 * t241 + (t232 * t255 + t234 * t249) * t230, t225, 0, 0; -t226, -t232 * t243 + (t231 * t249 - t232 * t253) * t230, -t227, 0, 0; 0, -t229 * t247 - t246, -t229 * t251 - t245, 0, 0; t227, t236 * t229 + t234 * t245, t226, 0, 0; t225, t237 * t229 + t231 * t245, t228, 0, 0; 0, t229 * t252 - t232 * t247, t229 * t248 - t232 * t251, 0, 0; -t237, -t231 * t254 - t242, 0, 0, 0; t236, t233 * t253 - t244, 0, 0, 0; 0, t251, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:28
	% EndTime: 2019-12-29 20:16:28
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (48->26), mult. (173->50), div. (0->0), fcn. (173->6), ass. (0->33)
	t297 = cos(qJ(1));
	t296 = cos(qJ(2));
	t316 = qJD(1) * t296;
	t301 = -qJD(3) + t316;
	t319 = t301 * t297;
	t294 = sin(qJ(1));
	t293 = sin(qJ(2));
	t312 = qJD(2) * t297;
	t304 = t293 * t312;
	t318 = t301 * t294 + t304;
	t317 = qJD(1) * t294;
	t315 = qJD(1) * t297;
	t314 = qJD(2) * t293;
	t313 = qJD(2) * t296;
	t292 = sin(qJ(3));
	t311 = qJD(3) * t292;
	t310 = qJD(3) * t293;
	t309 = qJD(3) * t296;
	t295 = cos(qJ(3));
	t308 = t295 * t314;
	t307 = t295 * t310;
	t306 = t294 * t314;
	t305 = t294 * t313;
	t303 = t296 * t312;
	t302 = qJD(1) - t309;
	t300 = t302 * t297;
	t299 = -t293 * t315 - t305;
	t298 = t293 * t317 - t303;
	t291 = t295 * t319 + (t302 * t292 - t308) * t294;
	t290 = t302 * t295 * t294 + (t306 - t319) * t292;
	t289 = t292 * t300 - t318 * t295;
	t288 = t318 * t292 + t295 * t300;
	t1 = [-t291, -t295 * t303 + (t295 * t317 + t297 * t311) * t293, t288, 0, 0; t289, -t295 * t305 + (t294 * t311 - t295 * t315) * t293, t290, 0, 0; 0, -t292 * t309 - t308, -t292 * t313 - t307, 0, 0; t299, -t294 * t316 - t304, 0, 0, 0; -t298, t296 * t315 - t306, 0, 0, 0; 0, t313, 0, 0, 0; t290, t298 * t292 - t297 * t307, t289, 0, 0; -t288, t299 * t292 - t294 * t307, t291, 0, 0; 0, -t292 * t314 + t295 * t309, -t292 * t310 + t295 * t313, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:16:29
	% EndTime: 2019-12-29 20:16:29
	% DurationCPUTime: 0.48s
	% Computational Cost: add. (219->57), mult. (709->92), div. (0->0), fcn. (733->8), ass. (0->47)
	t359 = sin(qJ(5));
	t360 = sin(qJ(3));
	t363 = cos(qJ(5));
	t364 = cos(qJ(3));
	t376 = t359 * t360 + t363 * t364;
	t401 = qJD(3) - qJD(5);
	t367 = t401 * t376;
	t377 = t359 * t364 - t360 * t363;
	t402 = t401 * t377;
	t365 = cos(qJ(2));
	t361 = sin(qJ(2));
	t362 = sin(qJ(1));
	t386 = t362 * qJD(2) * t361;
	t366 = cos(qJ(1));
	t391 = qJD(1) * t366;
	t400 = -t365 * t391 + t386;
	t394 = t366 * t364;
	t396 = t362 * t365;
	t345 = t360 * t396 + t394;
	t387 = qJD(3) * t366;
	t383 = t364 * t387;
	t389 = qJD(2) * t366;
	t384 = t361 * t389;
	t388 = qJD(3) * t362;
	t385 = t360 * t388;
	t341 = t345 * qJD(1) + t360 * t384 - t365 * t383 - t385;
	t392 = qJD(1) * t365;
	t395 = t366 * t360;
	t342 = (-qJD(3) * t365 + qJD(1)) * t395 + (-t384 + (qJD(3) - t392) * t362) * t364;
	t347 = -t362 * t364 + t365 * t395;
	t348 = t362 * t360 + t365 * t394;
	t399 = t341 * t363 + t342 * t359 + (t347 * t359 + t348 * t363) * qJD(5);
	t393 = qJD(1) * t362;
	t343 = (t388 * t365 - t393) * t364 + (-t387 - t400) * t360;
	t344 = t348 * qJD(1) - t364 * t386 - t365 * t385 - t383;
	t346 = t364 * t396 - t395;
	t398 = t343 * t359 + t344 * t363 + (t345 * t363 - t346 * t359) * qJD(5);
	t371 = (t345 * t359 + t346 * t363) * qJD(5) - t343 * t363 + t344 * t359;
	t390 = qJD(2) * t365;
	t375 = qJD(2) * t377;
	t374 = qJD(2) * t376;
	t373 = t365 * t374;
	t372 = t365 * t375;
	t337 = -t341 * t359 + t342 * t363 + (t347 * t363 - t348 * t359) * qJD(5);
	t339 = t361 * t402 + t376 * t390;
	t338 = -t367 * t361 + t372;
	t1 = [-t398, -t366 * t373 + (-t366 * t402 + t376 * t393) * t361, t399, 0, -t399; t337, -t362 * t373 + (-t362 * t402 - t376 * t391) * t361, t371, 0, -t371; 0, -t361 * t374 + t365 * t402, t338, 0, -t338; t371, t366 * t372 + (-t367 * t366 - t377 * t393) * t361, t337, 0, -t337; -t399, t362 * t372 + (-t367 * t362 + t377 * t391) * t361, t398, 0, -t398; 0, t361 * t375 + t365 * t367, t339, 0, -t339; t361 * t391 + t362 * t390, t362 * t392 + t384, 0, 0, 0; t361 * t393 - t365 * t389, t400, 0, 0, 0; 0, -t390, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end