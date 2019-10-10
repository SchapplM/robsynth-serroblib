% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRPR6_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR6_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (6->6), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = cos(pkin(6));
	t62 = sin(qJ(2));
	t66 = t61 * t62;
	t63 = cos(qJ(2));
	t65 = t61 * t63;
	t64 = qJD(2) * sin(pkin(6));
	t60 = cos(pkin(10));
	t58 = sin(pkin(10));
	t1 = [0, (t58 * t66 - t60 * t63) * qJD(2), 0, 0, 0, 0; 0, (-t58 * t63 - t60 * t66) * qJD(2), 0, 0, 0, 0; 0, -t62 * t64, 0, 0, 0, 0; 0, (t58 * t65 + t60 * t62) * qJD(2), 0, 0, 0, 0; 0, (t58 * t62 - t60 * t65) * qJD(2), 0, 0, 0, 0; 0, -t63 * t64, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:09
	% EndTime: 2019-10-09 21:39:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (4->4), mult. (24->17), div. (0->0), fcn. (24->6), ass. (0->9)
	t139 = cos(pkin(6));
	t140 = sin(qJ(2));
	t144 = t139 * t140;
	t141 = cos(qJ(2));
	t143 = t139 * t141;
	t142 = qJD(2) * sin(pkin(6));
	t138 = cos(pkin(10));
	t136 = sin(pkin(10));
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, (-t136 * t144 + t138 * t141) * qJD(2), 0, 0, 0, 0; 0, (t136 * t141 + t138 * t144) * qJD(2), 0, 0, 0, 0; 0, t140 * t142, 0, 0, 0, 0; 0, (-t136 * t143 - t138 * t140) * qJD(2), 0, 0, 0, 0; 0, (-t136 * t140 + t138 * t143) * qJD(2), 0, 0, 0, 0; 0, t141 * t142, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:10
	% EndTime: 2019-10-09 21:39:10
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (37->25), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t209 = sin(pkin(6));
	t212 = sin(qJ(4));
	t225 = t209 * t212;
	t214 = cos(qJ(4));
	t224 = t209 * t214;
	t215 = cos(qJ(2));
	t223 = t209 * t215;
	t211 = cos(pkin(6));
	t213 = sin(qJ(2));
	t222 = t211 * t213;
	t221 = t211 * t215;
	t220 = qJD(2) * t215;
	t219 = qJD(4) * t212;
	t218 = qJD(4) * t214;
	t217 = t209 * qJD(2) * t213;
	t208 = sin(pkin(10));
	t210 = cos(pkin(10));
	t216 = -t208 * t213 + t210 * t221;
	t205 = t208 * t215 + t210 * t222;
	t206 = t208 * t221 + t210 * t213;
	t207 = -t208 * t222 + t210 * t215;
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t201 = t205 * qJD(2);
	t200 = t216 * qJD(2);
	t1 = [0, -t202 * t212 + t207 * t218, 0, t203 * t214 + (-t206 * t212 - t208 * t224) * qJD(4), 0, 0; 0, t200 * t212 + t205 * t218, 0, t201 * t214 + (t210 * t224 + t212 * t216) * qJD(4), 0, 0; 0, (t212 * t220 + t213 * t218) * t209, 0, t214 * t217 + (-t211 * t214 + t212 * t223) * qJD(4), 0, 0; 0, -t202 * t214 - t207 * t219, 0, -t203 * t212 + (-t206 * t214 + t208 * t225) * qJD(4), 0, 0; 0, t200 * t214 - t205 * t219, 0, -t201 * t212 + (-t210 * t225 + t214 * t216) * qJD(4), 0, 0; 0, (-t213 * t219 + t214 * t220) * t209, 0, -t212 * t217 + (t211 * t212 + t214 * t223) * qJD(4), 0, 0; 0, -t203, 0, 0, 0, 0; 0, -t201, 0, 0, 0, 0; 0, -t217, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:11
	% EndTime: 2019-10-09 21:39:11
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (64->33), mult. (240->82), div. (0->0), fcn. (252->10), ass. (0->35)
	t304 = sin(pkin(6));
	t308 = sin(qJ(4));
	t325 = t304 * t308;
	t310 = cos(qJ(4));
	t324 = t304 * t310;
	t307 = cos(pkin(6));
	t309 = sin(qJ(2));
	t323 = t307 * t309;
	t311 = cos(qJ(2));
	t322 = t307 * t311;
	t321 = t308 * t311;
	t320 = qJD(2) * t311;
	t319 = qJD(4) * t308;
	t318 = qJD(4) * t310;
	t303 = sin(pkin(10));
	t317 = t303 * t323;
	t316 = qJD(2) * t304 * t309;
	t315 = t309 * t318;
	t306 = cos(pkin(10));
	t314 = -t303 * t309 + t306 * t322;
	t298 = t303 * t311 + t306 * t323;
	t299 = t303 * t322 + t306 * t309;
	t293 = t314 * qJD(2);
	t313 = t293 * t308 + t298 * t318;
	t295 = t299 * qJD(2);
	t300 = t306 * t311 - t317;
	t312 = -t295 * t308 + t300 * t318;
	t305 = cos(pkin(11));
	t302 = sin(pkin(11));
	t296 = -qJD(2) * t317 + t306 * t320;
	t294 = t298 * qJD(2);
	t292 = t310 * t316 + (t304 * t321 - t307 * t310) * qJD(4);
	t291 = t294 * t310 + (t306 * t324 + t308 * t314) * qJD(4);
	t290 = t296 * t310 + (-t299 * t308 - t303 * t324) * qJD(4);
	t1 = [0, -t296 * t302 + t312 * t305, 0, t290 * t305, 0, 0; 0, -t294 * t302 + t313 * t305, 0, t291 * t305, 0, 0; 0, (t305 * t315 + (-t302 * t309 + t305 * t321) * qJD(2)) * t304, 0, t292 * t305, 0, 0; 0, -t296 * t305 - t312 * t302, 0, -t290 * t302, 0, 0; 0, -t294 * t305 - t313 * t302, 0, -t291 * t302, 0, 0; 0, (-t302 * t315 + (-t302 * t321 - t305 * t309) * qJD(2)) * t304, 0, -t292 * t302, 0, 0; 0, t295 * t310 + t300 * t319, 0, t296 * t308 + (t299 * t310 - t303 * t325) * qJD(4), 0, 0; 0, -t293 * t310 + t298 * t319, 0, t294 * t308 + (t306 * t325 - t310 * t314) * qJD(4), 0, 0; 0, (t309 * t319 - t310 * t320) * t304, 0, t308 * t316 + (-t307 * t308 - t311 * t324) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:39:11
	% EndTime: 2019-10-09 21:39:11
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (219->60), mult. (516->126), div. (0->0), fcn. (564->10), ass. (0->53)
	t372 = sin(qJ(4));
	t373 = sin(qJ(2));
	t375 = cos(qJ(2));
	t374 = cos(qJ(4));
	t392 = qJD(4) * t374;
	t401 = (qJD(2) * t372 + qJD(6)) * t375 + t373 * t392;
	t369 = sin(pkin(6));
	t400 = t369 * t372;
	t399 = t369 * t373;
	t398 = t369 * t374;
	t397 = t369 * t375;
	t371 = cos(pkin(6));
	t396 = t371 * t373;
	t395 = t371 * t375;
	t394 = qJD(2) * t375;
	t393 = qJD(4) * t372;
	t367 = pkin(11) + qJ(6);
	t365 = sin(t367);
	t391 = qJD(6) * t365;
	t366 = cos(t367);
	t390 = qJD(6) * t366;
	t389 = qJD(6) * t372;
	t368 = sin(pkin(10));
	t388 = t368 * t396;
	t387 = qJD(2) * t399;
	t386 = t369 * t394;
	t370 = cos(pkin(10));
	t358 = t368 * t375 + t370 * t396;
	t354 = t358 * qJD(2);
	t383 = -t358 * t389 - t354;
	t356 = -qJD(2) * t388 + t370 * t394;
	t360 = t370 * t375 - t388;
	t382 = -t360 * t389 - t356;
	t381 = (-qJD(2) - t389) * t373;
	t379 = -t368 * t373 + t370 * t395;
	t380 = t370 * t398 + t372 * t379;
	t349 = t370 * t400 - t374 * t379;
	t359 = t368 * t395 + t370 * t373;
	t348 = t359 * t372 + t368 * t398;
	t347 = t359 * t374 - t368 * t400;
	t361 = -t371 * t372 - t374 * t397;
	t378 = -t371 * t374 + t372 * t397;
	t353 = t379 * qJD(2);
	t377 = qJD(6) * t379 + t353 * t372 + t358 * t392;
	t355 = t359 * qJD(2);
	t376 = -qJD(6) * t359 - t355 * t372 + t360 * t392;
	t352 = t378 * qJD(4) + t374 * t387;
	t351 = t361 * qJD(4) + t372 * t387;
	t346 = t380 * qJD(4) + t354 * t374;
	t345 = t349 * qJD(4) + t354 * t372;
	t344 = -t348 * qJD(4) + t356 * t374;
	t343 = t347 * qJD(4) + t356 * t372;
	t1 = [0, t382 * t365 + t376 * t366, 0, t344 * t366 - t347 * t391, 0, -t343 * t365 - t355 * t366 + (-t348 * t366 - t360 * t365) * qJD(6); 0, t383 * t365 + t377 * t366, 0, t346 * t366 - t349 * t391, 0, -t345 * t365 + t353 * t366 + (-t358 * t365 + t366 * t380) * qJD(6); 0, (t365 * t381 + t401 * t366) * t369, 0, t352 * t366 - t361 * t391, 0, t366 * t386 - t351 * t365 + (-t365 * t399 + t366 * t378) * qJD(6); 0, -t376 * t365 + t382 * t366, 0, -t344 * t365 - t347 * t390, 0, -t343 * t366 + t355 * t365 + (t348 * t365 - t360 * t366) * qJD(6); 0, -t377 * t365 + t383 * t366, 0, -t346 * t365 - t349 * t390, 0, -t345 * t366 - t353 * t365 + (-t358 * t366 - t365 * t380) * qJD(6); 0, (-t401 * t365 + t366 * t381) * t369, 0, -t352 * t365 - t361 * t390, 0, -t365 * t386 - t351 * t366 + (-t365 * t378 - t366 * t399) * qJD(6); 0, t355 * t374 + t360 * t393, 0, t343, 0, 0; 0, -t353 * t374 + t358 * t393, 0, t345, 0, 0; 0, (t373 * t393 - t374 * t394) * t369, 0, t351, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end