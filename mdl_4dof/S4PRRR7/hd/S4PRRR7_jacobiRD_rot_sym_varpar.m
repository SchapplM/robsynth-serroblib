% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% JRD_rot [9x4]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:35
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S4PRRR7_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_jacobiRD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_jacobiRD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR7_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:58
	% EndTime: 2019-12-29 12:34:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:34:59
	% EndTime: 2019-12-29 12:34:59
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(4));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(4));
	t60 = cos(pkin(8));
	t58 = sin(pkin(8));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0; 0, -t62 * t64, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0; 0, -t63 * t64, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:35:00
	% EndTime: 2019-12-29 12:35:00
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (36->22), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t210 = sin(pkin(4));
	t213 = sin(qJ(3));
	t226 = t210 * t213;
	t215 = cos(qJ(3));
	t225 = t210 * t215;
	t212 = cos(pkin(4));
	t214 = sin(qJ(2));
	t224 = t212 * t214;
	t216 = cos(qJ(2));
	t223 = t212 * t216;
	t222 = qJD(2) * t214;
	t221 = qJD(3) * t213;
	t220 = qJD(3) * t215;
	t219 = qJD(3) * t216;
	t218 = t210 * qJD(2) * t216;
	t209 = sin(pkin(8));
	t211 = cos(pkin(8));
	t205 = -t209 * t214 + t211 * t223;
	t206 = t209 * t216 + t211 * t224;
	t207 = -t209 * t223 - t211 * t214;
	t217 = t209 * t224 - t211 * t216;
	t204 = t217 * qJD(2);
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t201 = t205 * qJD(2);
	t1 = [0, t204 * t215 - t207 * t221, -t203 * t213 + (-t209 * t226 + t215 * t217) * qJD(3), 0; 0, -t202 * t215 - t205 * t221, -t201 * t213 + (-t206 * t215 + t211 * t226) * qJD(3), 0; 0, (-t213 * t219 - t215 * t222) * t210, -t213 * t218 + (-t212 * t213 - t214 * t225) * qJD(3), 0; 0, -t204 * t213 - t207 * t220, -t203 * t215 + (-t209 * t225 - t213 * t217) * qJD(3), 0; 0, t202 * t213 - t205 * t220, -t201 * t215 + (t206 * t213 + t211 * t225) * qJD(3), 0; 0, (t213 * t222 - t215 * t219) * t210, -t215 * t218 + (-t212 * t215 + t214 * t226) * qJD(3), 0; 0, t203, 0, 0; 0, t201, 0, 0; 0, t218, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:35:02
	% EndTime: 2019-12-29 12:35:02
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (153->59), mult. (516->127), div. (0->0), fcn. (564->10), ass. (0->53)
	t353 = sin(qJ(3));
	t354 = sin(qJ(2));
	t356 = cos(qJ(3));
	t357 = cos(qJ(2));
	t373 = qJD(3) * t357;
	t383 = (qJD(2) * t356 - qJD(4)) * t354 + t353 * t373;
	t349 = sin(pkin(4));
	t382 = t349 * t353;
	t381 = t349 * t356;
	t380 = t349 * t357;
	t351 = cos(pkin(4));
	t379 = t351 * t354;
	t378 = t351 * t357;
	t377 = qJD(2) * t354;
	t376 = qJD(2) * t357;
	t375 = qJD(3) * t353;
	t374 = qJD(3) * t356;
	t352 = sin(qJ(4));
	t372 = qJD(4) * t352;
	t355 = cos(qJ(4));
	t371 = qJD(4) * t355;
	t370 = qJD(4) * t356;
	t348 = sin(pkin(8));
	t369 = t348 * t379;
	t368 = t349 * t377;
	t367 = t349 * t376;
	t350 = cos(pkin(8));
	t360 = -t348 * t354 + t350 * t378;
	t336 = t360 * qJD(2);
	t364 = -t360 * t370 + t336;
	t342 = t348 * t378 + t350 * t354;
	t338 = t342 * qJD(2);
	t363 = t342 * t370 - t338;
	t362 = (qJD(2) - t370) * t357;
	t341 = t348 * t357 + t350 * t379;
	t330 = -t341 * t353 - t350 * t381;
	t361 = -t341 * t356 + t350 * t382;
	t343 = t350 * t357 - t369;
	t332 = -t343 * t353 + t348 * t381;
	t333 = t343 * t356 + t348 * t382;
	t345 = t351 * t353 + t354 * t381;
	t344 = t351 * t356 - t354 * t382;
	t337 = t341 * qJD(2);
	t359 = qJD(4) * t341 - t337 * t356 - t360 * t375;
	t339 = -qJD(2) * t369 + t350 * t376;
	t358 = qJD(4) * t343 - t339 * t356 + t342 * t375;
	t335 = t344 * qJD(3) + t356 * t367;
	t334 = -t345 * qJD(3) - t353 * t367;
	t329 = t332 * qJD(3) - t338 * t356;
	t328 = -t333 * qJD(3) + t338 * t353;
	t327 = t330 * qJD(3) + t336 * t356;
	t326 = t361 * qJD(3) - t336 * t353;
	t1 = [0, t363 * t352 + t358 * t355, t328 * t355 - t332 * t372, -t329 * t352 + t339 * t355 + (-t333 * t355 - t342 * t352) * qJD(4); 0, t364 * t352 + t359 * t355, t326 * t355 - t330 * t372, -t327 * t352 + t337 * t355 + (t352 * t360 + t355 * t361) * qJD(4); 0, (t352 * t362 - t383 * t355) * t349, t334 * t355 - t344 * t372, t355 * t368 - t335 * t352 + (-t345 * t355 + t352 * t380) * qJD(4); 0, -t358 * t352 + t363 * t355, -t328 * t352 - t332 * t371, -t329 * t355 - t339 * t352 + (t333 * t352 - t342 * t355) * qJD(4); 0, -t359 * t352 + t364 * t355, -t326 * t352 - t330 * t371, -t327 * t355 - t337 * t352 + (-t352 * t361 + t355 * t360) * qJD(4); 0, (t383 * t352 + t355 * t362) * t349, -t334 * t352 - t344 * t371, -t352 * t368 - t335 * t355 + (t345 * t352 + t355 * t380) * qJD(4); 0, -t339 * t353 - t342 * t374, t329, 0; 0, -t337 * t353 + t360 * t374, t327, 0; 0, (-t353 * t377 + t356 * t373) * t349, t335, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,4);
end