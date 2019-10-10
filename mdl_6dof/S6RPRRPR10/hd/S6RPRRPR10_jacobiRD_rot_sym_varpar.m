% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RPRRPR10_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR10_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:39:30
	% EndTime: 2019-10-10 01:39:30
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:39:30
	% EndTime: 2019-10-10 01:39:30
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
	% StartTime: 2019-10-10 01:39:30
	% EndTime: 2019-10-10 01:39:30
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t13 = qJD(1) * sin(qJ(1));
	t11 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; t13, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t13, 0, 0, 0, 0, 0; t11, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:39:30
	% EndTime: 2019-10-10 01:39:30
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
	% StartTime: 2019-10-10 01:39:31
	% EndTime: 2019-10-10 01:39:31
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
	% StartTime: 2019-10-10 01:39:32
	% EndTime: 2019-10-10 01:39:32
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (49->25), mult. (173->45), div. (0->0), fcn. (173->6), ass. (0->32)
	t280 = sin(qJ(1));
	t279 = sin(qJ(3));
	t289 = qJD(1) * t279 + qJD(4);
	t282 = cos(qJ(3));
	t283 = cos(qJ(1));
	t299 = qJD(3) * t283;
	t291 = t282 * t299;
	t304 = t289 * t280 - t291;
	t303 = qJD(1) * t280;
	t302 = qJD(1) * t283;
	t301 = qJD(3) * t279;
	t300 = qJD(3) * t282;
	t298 = qJD(4) * t279;
	t297 = qJD(4) * t282;
	t296 = t282 * t302;
	t281 = cos(qJ(4));
	t295 = t281 * t300;
	t294 = t280 * t300;
	t278 = sin(qJ(4));
	t293 = t278 * t297;
	t292 = t281 * t297;
	t290 = qJD(1) + t298;
	t288 = t290 * t283;
	t287 = t289 * t283;
	t286 = -t280 * t301 + t296;
	t285 = t279 * t299 + t282 * t303;
	t284 = -t281 * t301 - t293;
	t277 = t281 * t287 + (-t290 * t278 + t295) * t280;
	t276 = t290 * t281 * t280 + (t287 + t294) * t278;
	t275 = t278 * t288 + t304 * t281;
	t274 = -t304 * t278 + t281 * t288;
	t1 = [-t275, 0, t284 * t280 + t281 * t296, -t276, 0, 0; t277, 0, t285 * t281 + t283 * t293, t274, 0, 0; 0, 0, t278 * t298 - t295, t278 * t301 - t292, 0, 0; t285, 0, t279 * t302 + t294, 0, 0, 0; -t286, 0, t279 * t303 - t291, 0, 0, 0; 0, 0, -t301, 0, 0, 0; t274, 0, t286 * t278 + t280 * t292, t277, 0, 0; t276, 0, t285 * t278 - t283 * t292, t275, 0, 0; 0, 0, -t278 * t300 - t281 * t298, t284, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:39:32
	% EndTime: 2019-10-10 01:39:33
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (218->54), mult. (709->93), div. (0->0), fcn. (733->8), ass. (0->51)
	t399 = qJD(4) - qJD(6);
	t355 = sin(qJ(4));
	t359 = cos(qJ(4));
	t356 = sin(qJ(3));
	t380 = qJD(1) * t356 + qJD(4);
	t360 = cos(qJ(3));
	t361 = cos(qJ(1));
	t385 = qJD(3) * t361;
	t382 = t360 * t385;
	t390 = t361 * t359;
	t384 = t356 * t390;
	t388 = qJD(1) * t361;
	t357 = sin(qJ(1));
	t394 = t357 * t355;
	t340 = -qJD(4) * t384 - t355 * t382 - t359 * t388 + t380 * t394;
	t391 = t361 * t355;
	t393 = t357 * t359;
	t345 = t356 * t393 + t391;
	t346 = t356 * t391 + t393;
	t341 = t345 * qJD(1) + t346 * qJD(4) - t359 * t382;
	t347 = t384 - t394;
	t354 = sin(qJ(6));
	t358 = cos(qJ(6));
	t364 = -t340 * t354 - t341 * t358 + (t346 * t358 - t347 * t354) * qJD(6);
	t372 = t380 * t361;
	t381 = qJD(4) * t356 + qJD(1);
	t386 = qJD(3) * t360;
	t383 = t357 * t386;
	t342 = t381 * t393 + (t372 + t383) * t355;
	t343 = t359 * t372 + (-t381 * t355 + t359 * t386) * t357;
	t344 = t356 * t394 - t390;
	t398 = -t342 * t358 + t343 * t354 + (t344 * t354 + t345 * t358) * qJD(6);
	t397 = t399 * t360;
	t396 = t340 * t358 - t341 * t354 + (t346 * t354 + t347 * t358) * qJD(6);
	t392 = t360 * t361;
	t389 = qJD(1) * t357;
	t387 = qJD(3) * t356;
	t374 = t354 * t359 - t355 * t358;
	t373 = t354 * t355 + t358 * t359;
	t371 = qJD(1) * t374;
	t370 = qJD(1) * t373;
	t369 = qJD(3) * t374;
	t368 = qJD(3) * t373;
	t367 = t374 * t387;
	t366 = t373 * t387;
	t336 = t342 * t354 + t343 * t358 + (t344 * t358 - t345 * t354) * qJD(6);
	t363 = t399 * t374;
	t362 = t399 * t373;
	t338 = -t363 * t360 + t366;
	t337 = -t362 * t360 - t367;
	t1 = [t364, 0, t370 * t392 + (-t356 * t368 + t374 * t397) * t357, t398, 0, -t398; t336, 0, t361 * t366 + (t357 * t370 - t363 * t361) * t360, -t396, 0, t396; 0, 0, -t363 * t356 - t360 * t368, t337, 0, -t337; -t396, 0, -t371 * t392 + (t356 * t369 + t373 * t397) * t357, t336, 0, -t336; -t398, 0, -t361 * t367 + (-t357 * t371 - t362 * t361) * t360, -t364, 0, t364; 0, 0, -t362 * t356 + t360 * t369, -t338, 0, t338; -t356 * t385 - t360 * t389, 0, -t356 * t388 - t383, 0, 0, 0; -t357 * t387 + t360 * t388, 0, -t356 * t389 + t382, 0, 0, 0; 0, 0, t387, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end