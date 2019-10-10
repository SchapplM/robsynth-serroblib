% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:48
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6RRPPRR8_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR8_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobiRD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
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
	% StartTime: 2019-10-10 09:48:31
	% EndTime: 2019-10-10 09:48:31
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
	% StartTime: 2019-10-10 09:48:32
	% EndTime: 2019-10-10 09:48:33
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
	% StartTime: 2019-10-10 09:48:32
	% EndTime: 2019-10-10 09:48:32
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (18->15), mult. (77->37), div. (0->0), fcn. (77->6), ass. (0->18)
	t215 = sin(qJ(1));
	t216 = cos(qJ(2));
	t228 = t215 * t216;
	t217 = cos(qJ(1));
	t227 = t216 * t217;
	t226 = qJD(1) * t215;
	t225 = qJD(1) * t217;
	t214 = sin(qJ(2));
	t224 = qJD(2) * t214;
	t223 = qJD(2) * t216;
	t222 = qJD(2) * t217;
	t221 = t215 * t224;
	t220 = t214 * t222;
	t219 = -t214 * t225 - t215 * t223;
	t218 = t214 * t226 - t216 * t222;
	t213 = cos(pkin(10));
	t212 = sin(pkin(10));
	t1 = [t213 * t221 + (-t212 * t215 - t213 * t227) * qJD(1), t218 * t213, 0, 0, 0, 0; -t213 * t220 + (t212 * t217 - t213 * t228) * qJD(1), t219 * t213, 0, 0, 0, 0; 0, -t213 * t224, 0, 0, 0, 0; t219, -t216 * t226 - t220, 0, 0, 0, 0; -t218, t216 * t225 - t221, 0, 0, 0, 0; 0, t223, 0, 0, 0, 0; t212 * t221 + (-t212 * t227 + t213 * t215) * qJD(1), t218 * t212, 0, 0, 0, 0; -t212 * t220 + (-t212 * t228 - t213 * t217) * qJD(1), t219 * t212, 0, 0, 0, 0; 0, -t212 * t224, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:33
	% EndTime: 2019-10-10 09:48:33
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (109->41), mult. (385->83), div. (0->0), fcn. (401->8), ass. (0->40)
	t304 = cos(qJ(2));
	t301 = sin(qJ(2));
	t302 = sin(qJ(1));
	t320 = t302 * qJD(2) * t301;
	t305 = cos(qJ(1));
	t324 = qJD(1) * t305;
	t331 = -t304 * t324 + t320;
	t298 = sin(pkin(10));
	t299 = cos(pkin(10));
	t300 = sin(qJ(5));
	t303 = cos(qJ(5));
	t315 = t298 * t300 + t299 * t303;
	t330 = qJD(5) * t315;
	t325 = qJD(1) * t302;
	t288 = -t331 * t298 - t299 * t325;
	t326 = t305 * t299;
	t293 = t302 * t298 + t304 * t326;
	t289 = t293 * qJD(1) - t299 * t320;
	t328 = t302 * t304;
	t290 = t298 * t328 + t326;
	t327 = t305 * t298;
	t291 = t299 * t328 - t327;
	t329 = -t288 * t303 + t289 * t300 + (t290 * t300 + t291 * t303) * qJD(5);
	t323 = qJD(2) * t304;
	t322 = qJD(2) * t305;
	t319 = t301 * t322;
	t316 = t298 * t303 - t299 * t300;
	t314 = t316 * t304;
	t313 = qJD(2) * t316;
	t312 = qJD(2) * t315;
	t311 = qJD(5) * t316;
	t308 = t304 * t313;
	t307 = t304 * t312;
	t306 = -t288 * t300 - t289 * t303 + (-t290 * t303 + t291 * t300) * qJD(5);
	t292 = -t302 * t299 + t304 * t327;
	t287 = -t291 * qJD(1) - t299 * t319;
	t286 = -t290 * qJD(1) - t298 * t319;
	t285 = t286 * t300 + t287 * t303 + (t292 * t303 - t293 * t300) * qJD(5);
	t284 = t286 * t303 - t287 * t300 + (-t292 * t300 - t293 * t303) * qJD(5);
	t1 = [t306, -t305 * t307 + (-t305 * t311 + t315 * t325) * t301, 0, 0, t284, 0; t285, -t302 * t307 + (-t302 * t311 - t315 * t324) * t301, 0, 0, -t329, 0; 0, qJD(5) * t314 - t301 * t312, 0, 0, qJD(2) * t314 - t301 * t330, 0; t329, -t305 * t308 + (t305 * t330 + t316 * t325) * t301, 0, 0, -t285, 0; t284, -t302 * t308 + (t302 * t330 - t316 * t324) * t301, 0, 0, t306, 0; 0, -t301 * t313 - t304 * t330, 0, 0, -t301 * t311 - t307, 0; t301 * t324 + t302 * t323, t304 * t325 + t319, 0, 0, 0, 0; t301 * t325 - t304 * t322, t331, 0, 0, 0, 0; 0, -t323, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:48:34
	% EndTime: 2019-10-10 09:48:34
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (291->41), mult. (521->77), div. (0->0), fcn. (545->8), ass. (0->45)
	t358 = qJD(5) + qJD(6);
	t359 = qJ(5) + qJ(6);
	t356 = sin(t359);
	t357 = cos(t359);
	t360 = sin(pkin(10));
	t361 = cos(pkin(10));
	t377 = t356 * t360 + t357 * t361;
	t374 = t358 * t377;
	t365 = cos(qJ(1));
	t389 = t365 * t361;
	t363 = sin(qJ(1));
	t364 = cos(qJ(2));
	t391 = t363 * t364;
	t349 = t360 * t391 + t389;
	t376 = t363 * t360 + t364 * t389;
	t362 = sin(qJ(2));
	t385 = qJD(2) * t365;
	t382 = t362 * t385;
	t394 = -t349 * qJD(1) - t376 * t358 - t360 * t382;
	t383 = t363 * qJD(2) * t362;
	t387 = qJD(1) * t365;
	t393 = -t364 * t387 + t383;
	t379 = t376 * qJD(1) + t349 * t358 - t361 * t383;
	t390 = t365 * t360;
	t350 = t361 * t391 - t390;
	t388 = qJD(1) * t363;
	t380 = t350 * t358 + t393 * t360 + t361 * t388;
	t392 = t379 * t356 + t380 * t357;
	t386 = qJD(2) * t364;
	t381 = -t350 * qJD(1) + (-t363 * t361 + t364 * t390) * t358 - t361 * t382;
	t378 = t356 * t361 - t357 * t360;
	t375 = t358 * t378;
	t373 = t378 * t364;
	t372 = t378 * t365;
	t371 = t377 * t363;
	t370 = qJD(2) * t378;
	t369 = qJD(2) * t377;
	t367 = t364 * t370;
	t366 = t364 * t369;
	t341 = t380 * t356 - t379 * t357;
	t343 = t362 * t375 - t366;
	t342 = -qJD(2) * t373 - t362 * t374;
	t339 = t394 * t356 + t381 * t357;
	t338 = -t381 * t356 + t394 * t357;
	t1 = [t341, -t365 * t366 + (qJD(1) * t371 + t358 * t372) * t362, 0, 0, t338, t338; t339, -t363 * t366 + (t363 * t375 - t377 * t387) * t362, 0, 0, -t392, -t392; 0, -t358 * t373 - t362 * t369, 0, 0, t342, t342; t392, t365 * t367 + (t365 * t374 - t378 * t388) * t362, 0, 0, -t339, -t339; t338, t363 * t367 + (qJD(1) * t372 + t358 * t371) * t362, 0, 0, t341, t341; 0, t362 * t370 - t364 * t374, 0, 0, t343, t343; t362 * t387 + t363 * t386, t364 * t388 + t382, 0, 0, 0, 0; t362 * t388 - t364 * t385, t393, 0, 0, 0, 0; 0, -t386, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end