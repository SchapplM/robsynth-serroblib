% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:17
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPRPR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
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
	% StartTime: 2019-10-10 10:17:10
	% EndTime: 2019-10-10 10:17:10
	% DurationCPUTime: 0.04s
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
	t1 = [t32, t29, 0, 0, 0, 0; -t30, -t31, 0, 0, 0, 0; 0, -t39, 0, 0, 0, 0; t31, t30, 0, 0, 0, 0; t29, t32, 0, 0, 0, 0; 0, -t38, 0, 0, 0, 0; -t41, 0, 0, 0, 0, 0; t40, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:11
	% EndTime: 2019-10-10 10:17:11
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (17->17), mult. (77->37), div. (0->0), fcn. (77->6), ass. (0->18)
	t180 = sin(qJ(1));
	t181 = cos(qJ(2));
	t193 = t180 * t181;
	t182 = cos(qJ(1));
	t192 = t181 * t182;
	t191 = qJD(1) * t180;
	t190 = qJD(1) * t182;
	t179 = sin(qJ(2));
	t189 = qJD(2) * t179;
	t188 = qJD(2) * t181;
	t187 = qJD(2) * t182;
	t186 = t180 * t189;
	t185 = t179 * t187;
	t184 = t179 * t190 + t180 * t188;
	t183 = t179 * t191 - t181 * t187;
	t178 = cos(pkin(10));
	t177 = sin(pkin(10));
	t1 = [t178 * t186 + (-t177 * t180 - t178 * t192) * qJD(1), t183 * t178, 0, 0, 0, 0; -t178 * t185 + (t177 * t182 - t178 * t193) * qJD(1), -t184 * t178, 0, 0, 0, 0; 0, -t178 * t189, 0, 0, 0, 0; -t177 * t186 + (t177 * t192 - t178 * t180) * qJD(1), -t183 * t177, 0, 0, 0, 0; t177 * t185 + (t177 * t193 + t178 * t182) * qJD(1), t184 * t177, 0, 0, 0, 0; 0, t177 * t189, 0, 0, 0, 0; -t184, -t181 * t191 - t185, 0, 0, 0, 0; -t183, t181 * t190 - t186, 0, 0, 0, 0; 0, t188, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:11
	% EndTime: 2019-10-10 10:17:12
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (108->25), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->31)
	t253 = cos(qJ(1));
	t252 = cos(qJ(2));
	t264 = qJD(4) * t252;
	t259 = -qJD(1) + t264;
	t272 = t253 * t259;
	t270 = qJD(1) * t252;
	t258 = -qJD(4) + t270;
	t251 = sin(qJ(1));
	t250 = sin(qJ(2));
	t268 = qJD(2) * t250;
	t261 = t251 * t268;
	t271 = t258 * t253 - t261;
	t269 = qJD(1) * t253;
	t267 = qJD(2) * t252;
	t266 = qJD(2) * t253;
	t265 = qJD(4) * t250;
	t249 = pkin(10) + qJ(4);
	t247 = sin(t249);
	t263 = t247 * t265;
	t248 = cos(t249);
	t262 = t248 * t265;
	t260 = t250 * t266;
	t257 = t259 * t251;
	t256 = t250 * t269 + t251 * t267;
	t255 = qJD(1) * t251 * t250 - t252 * t266;
	t254 = t258 * t251 + t260;
	t246 = t247 * t257 - t271 * t248;
	t245 = t271 * t247 + t248 * t257;
	t244 = t247 * t272 + t254 * t248;
	t243 = t254 * t247 - t248 * t272;
	t1 = [t246, t255 * t248 + t253 * t263, 0, t243, 0, 0; -t244, -t256 * t248 + t251 * t263, 0, -t245, 0, 0; 0, -t247 * t264 - t248 * t268, 0, -t247 * t267 - t262, 0, 0; t245, -t255 * t247 + t253 * t262, 0, t244, 0, 0; t243, t256 * t247 + t251 * t262, 0, t246, 0, 0; 0, t247 * t268 - t248 * t264, 0, -t248 * t267 + t263, 0, 0; -t256, -t251 * t270 - t260, 0, 0, 0, 0; -t255, t252 * t269 - t261, 0, 0, 0, 0; 0, t267, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:12
	% EndTime: 2019-10-10 10:17:12
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (108->24), mult. (173->44), div. (0->0), fcn. (173->6), ass. (0->31)
	t311 = sin(qJ(1));
	t312 = cos(qJ(2));
	t330 = qJD(1) * t312;
	t318 = -qJD(4) + t330;
	t310 = sin(qJ(2));
	t313 = cos(qJ(1));
	t326 = qJD(2) * t313;
	t320 = t310 * t326;
	t332 = t318 * t311 + t320;
	t328 = qJD(2) * t310;
	t321 = t311 * t328;
	t331 = t318 * t313 - t321;
	t329 = qJD(1) * t313;
	t327 = qJD(2) * t312;
	t325 = qJD(4) * t310;
	t324 = qJD(4) * t312;
	t309 = pkin(10) + qJ(4);
	t307 = sin(t309);
	t323 = t307 * t325;
	t308 = cos(t309);
	t322 = t308 * t325;
	t319 = qJD(1) - t324;
	t317 = t319 * t311;
	t316 = t319 * t313;
	t315 = -t310 * t329 - t311 * t327;
	t314 = qJD(1) * t311 * t310 - t312 * t326;
	t306 = t307 * t317 + t331 * t308;
	t305 = -t331 * t307 + t308 * t317;
	t304 = t307 * t316 - t332 * t308;
	t303 = t332 * t307 + t308 * t316;
	t1 = [-t306, t314 * t308 + t313 * t323, 0, t303, 0, 0; t304, t315 * t308 + t311 * t323, 0, t305, 0, 0; 0, -t307 * t324 - t308 * t328, 0, -t307 * t327 - t322, 0, 0; t315, -t311 * t330 - t320, 0, 0, 0, 0; -t314, t312 * t329 - t321, 0, 0, 0, 0; 0, t327, 0, 0, 0, 0; t305, t314 * t307 - t313 * t322, 0, t304, 0, 0; -t303, t315 * t307 - t311 * t322, 0, t306, 0, 0; 0, -t307 * t328 + t308 * t324, 0, t308 * t327 - t323, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:17:13
	% EndTime: 2019-10-10 10:17:13
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (455->58), mult. (709->92), div. (0->0), fcn. (733->8), ass. (0->48)
	t384 = pkin(10) + qJ(4);
	t382 = sin(t384);
	t383 = cos(t384);
	t385 = sin(qJ(6));
	t388 = cos(qJ(6));
	t400 = t382 * t385 + t383 * t388;
	t425 = qJD(4) - qJD(6);
	t391 = t425 * t400;
	t401 = t382 * t388 - t383 * t385;
	t426 = t425 * t401;
	t389 = cos(qJ(2));
	t386 = sin(qJ(2));
	t387 = sin(qJ(1));
	t410 = t387 * qJD(2) * t386;
	t390 = cos(qJ(1));
	t415 = qJD(1) * t390;
	t424 = -t389 * t415 + t410;
	t418 = t390 * t383;
	t420 = t387 * t389;
	t368 = t382 * t420 + t418;
	t411 = qJD(4) * t390;
	t407 = t383 * t411;
	t413 = qJD(2) * t390;
	t408 = t386 * t413;
	t412 = qJD(4) * t387;
	t409 = t382 * t412;
	t364 = t368 * qJD(1) + t382 * t408 - t389 * t407 - t409;
	t416 = qJD(1) * t389;
	t419 = t390 * t382;
	t365 = (-qJD(4) * t389 + qJD(1)) * t419 + (-t408 + (qJD(4) - t416) * t387) * t383;
	t370 = -t387 * t383 + t389 * t419;
	t371 = t387 * t382 + t389 * t418;
	t423 = t364 * t388 + t365 * t385 + (t370 * t385 + t371 * t388) * qJD(6);
	t417 = qJD(1) * t387;
	t366 = (t412 * t389 - t417) * t383 + (-t411 - t424) * t382;
	t367 = t371 * qJD(1) - t383 * t410 - t389 * t409 - t407;
	t369 = t383 * t420 - t419;
	t422 = t366 * t385 + t367 * t388 + (t368 * t388 - t369 * t385) * qJD(6);
	t395 = (t368 * t385 + t369 * t388) * qJD(6) - t366 * t388 + t367 * t385;
	t414 = qJD(2) * t389;
	t399 = qJD(2) * t401;
	t398 = qJD(2) * t400;
	t397 = t389 * t398;
	t396 = t389 * t399;
	t360 = -t364 * t385 + t365 * t388 + (t370 * t388 - t371 * t385) * qJD(6);
	t362 = -t386 * t426 + t400 * t414;
	t361 = -t391 * t386 - t396;
	t1 = [-t422, -t390 * t397 + (t390 * t426 + t400 * t417) * t386, 0, t423, 0, -t423; t360, -t387 * t397 + (t387 * t426 - t400 * t415) * t386, 0, t395, 0, -t395; 0, -t386 * t398 - t389 * t426, 0, t361, 0, -t361; t395, -t390 * t396 + (-t391 * t390 + t401 * t417) * t386, 0, t360, 0, -t360; -t423, -t387 * t396 + (-t391 * t387 - t401 * t415) * t386, 0, t422, 0, -t422; 0, -t386 * t399 + t389 * t391, 0, t362, 0, -t362; t386 * t415 + t387 * t414, t387 * t416 + t408, 0, 0, 0, 0; t386 * t417 - t389 * t413, t424, 0, 0, 0, 0; 0, -t414, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end