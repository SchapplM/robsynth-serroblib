% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRPPR1_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR1_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:15
	% EndTime: 2019-10-09 22:07:15
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-09 22:07:16
	% EndTime: 2019-10-09 22:07:16
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (36->22), mult. (140->59), div. (0->0), fcn. (148->8), ass. (0->26)
	t211 = sin(pkin(6));
	t214 = sin(qJ(3));
	t227 = t211 * t214;
	t216 = cos(qJ(3));
	t226 = t211 * t216;
	t213 = cos(pkin(6));
	t215 = sin(qJ(2));
	t225 = t213 * t215;
	t217 = cos(qJ(2));
	t224 = t213 * t217;
	t223 = qJD(2) * t215;
	t222 = qJD(3) * t214;
	t221 = qJD(3) * t216;
	t220 = qJD(3) * t217;
	t219 = t211 * qJD(2) * t217;
	t210 = sin(pkin(10));
	t212 = cos(pkin(10));
	t206 = -t210 * t215 + t212 * t224;
	t207 = t210 * t217 + t212 * t225;
	t208 = -t210 * t224 - t212 * t215;
	t218 = t210 * t225 - t212 * t217;
	t205 = t218 * qJD(2);
	t204 = t208 * qJD(2);
	t203 = t207 * qJD(2);
	t202 = t206 * qJD(2);
	t1 = [0, t205 * t216 - t208 * t222, -t204 * t214 + (-t210 * t227 + t216 * t218) * qJD(3), 0, 0, 0; 0, -t203 * t216 - t206 * t222, -t202 * t214 + (-t207 * t216 + t212 * t227) * qJD(3), 0, 0, 0; 0, (-t214 * t220 - t216 * t223) * t211, -t214 * t219 + (-t213 * t214 - t215 * t226) * qJD(3), 0, 0, 0; 0, -t205 * t214 - t208 * t221, -t204 * t216 + (-t210 * t226 - t214 * t218) * qJD(3), 0, 0, 0; 0, t203 * t214 - t206 * t221, -t202 * t216 + (t207 * t214 + t212 * t226) * qJD(3), 0, 0, 0; 0, (t214 * t223 - t216 * t220) * t211, -t216 * t219 + (-t213 * t216 + t215 * t227) * qJD(3), 0, 0, 0; 0, t204, 0, 0, 0, 0; 0, t202, 0, 0, 0, 0; 0, t219, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:16
	% EndTime: 2019-10-09 22:07:16
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (66->23), mult. (140->60), div. (0->0), fcn. (148->8), ass. (0->28)
	t230 = sin(pkin(10));
	t231 = sin(pkin(6));
	t246 = t230 * t231;
	t232 = cos(pkin(10));
	t245 = t231 * t232;
	t234 = sin(qJ(2));
	t244 = t231 * t234;
	t233 = cos(pkin(6));
	t243 = t233 * t234;
	t235 = cos(qJ(2));
	t242 = t233 * t235;
	t241 = qJD(2) * t234;
	t229 = qJ(3) + pkin(11);
	t227 = sin(t229);
	t240 = qJD(3) * t227;
	t228 = cos(t229);
	t239 = qJD(3) * t228;
	t238 = qJD(3) * t235;
	t237 = t231 * qJD(2) * t235;
	t223 = -t230 * t234 + t232 * t242;
	t224 = t230 * t235 + t232 * t243;
	t225 = -t230 * t242 - t232 * t234;
	t236 = t230 * t243 - t232 * t235;
	t222 = t236 * qJD(2);
	t221 = t225 * qJD(2);
	t220 = t224 * qJD(2);
	t219 = t223 * qJD(2);
	t1 = [0, t222 * t228 - t225 * t240, -t221 * t227 + (-t227 * t246 + t228 * t236) * qJD(3), 0, 0, 0; 0, -t220 * t228 - t223 * t240, -t219 * t227 + (-t224 * t228 + t227 * t245) * qJD(3), 0, 0, 0; 0, (-t227 * t238 - t228 * t241) * t231, -t227 * t237 + (-t227 * t233 - t228 * t244) * qJD(3), 0, 0, 0; 0, -t222 * t227 - t225 * t239, -t221 * t228 + (-t227 * t236 - t228 * t246) * qJD(3), 0, 0, 0; 0, t220 * t227 - t223 * t239, -t219 * t228 + (t224 * t227 + t228 * t245) * qJD(3), 0, 0, 0; 0, (t227 * t241 - t228 * t238) * t231, -t228 * t237 + (t227 * t244 - t228 * t233) * qJD(3), 0, 0, 0; 0, t221, 0, 0, 0, 0; 0, t219, 0, 0, 0, 0; 0, t237, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:17
	% EndTime: 2019-10-09 22:07:17
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (109->33), mult. (240->83), div. (0->0), fcn. (252->10), ass. (0->36)
	t322 = qJ(3) + pkin(11);
	t321 = cos(t322);
	t329 = sin(qJ(2));
	t344 = t321 * t329;
	t324 = sin(pkin(10));
	t325 = sin(pkin(6));
	t343 = t324 * t325;
	t327 = cos(pkin(10));
	t342 = t325 * t327;
	t341 = t325 * t329;
	t328 = cos(pkin(6));
	t340 = t328 * t329;
	t330 = cos(qJ(2));
	t339 = t328 * t330;
	t320 = sin(t322);
	t338 = qJD(3) * t320;
	t337 = qJD(3) * t321;
	t336 = qJD(3) * t330;
	t335 = qJD(2) * t325 * t330;
	t334 = t320 * t336;
	t316 = -t324 * t329 + t327 * t339;
	t317 = t324 * t330 + t327 * t340;
	t318 = -t324 * t339 - t327 * t329;
	t333 = t324 * t340 - t327 * t330;
	t313 = t317 * qJD(2);
	t332 = t313 * t321 + t316 * t338;
	t315 = t333 * qJD(2);
	t331 = -t315 * t321 + t318 * t338;
	t326 = cos(pkin(12));
	t323 = sin(pkin(12));
	t314 = t318 * qJD(2);
	t312 = t316 * qJD(2);
	t311 = -t320 * t335 + (-t320 * t328 - t321 * t341) * qJD(3);
	t310 = -t314 * t320 + (-t320 * t343 + t321 * t333) * qJD(3);
	t309 = -t312 * t320 + (-t317 * t321 + t320 * t342) * qJD(3);
	t1 = [0, t314 * t323 - t331 * t326, t310 * t326, 0, 0, 0; 0, t312 * t323 - t332 * t326, t309 * t326, 0, 0, 0; 0, (-t326 * t334 + (t323 * t330 - t326 * t344) * qJD(2)) * t325, t311 * t326, 0, 0, 0; 0, t314 * t326 + t331 * t323, -t310 * t323, 0, 0, 0; 0, t312 * t326 + t332 * t323, -t309 * t323, 0, 0, 0; 0, (t323 * t334 + (t323 * t344 + t326 * t330) * qJD(2)) * t325, -t311 * t323, 0, 0, 0; 0, t315 * t320 + t318 * t337, t314 * t321 + (t320 * t333 + t321 * t343) * qJD(3), 0, 0, 0; 0, -t313 * t320 + t316 * t337, t312 * t321 + (-t317 * t320 - t321 * t342) * qJD(3), 0, 0, 0; 0, (-qJD(2) * t320 * t329 + t321 * t336) * t325, t321 * t335 + (-t320 * t341 + t321 * t328) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:07:18
	% EndTime: 2019-10-09 22:07:18
	% DurationCPUTime: 0.40s
	% Computational Cost: add. (312->61), mult. (516->128), div. (0->0), fcn. (564->10), ass. (0->56)
	t387 = qJ(3) + pkin(11);
	t383 = sin(t387);
	t385 = cos(t387);
	t392 = sin(qJ(2));
	t393 = cos(qJ(2));
	t409 = qJD(3) * t393;
	t420 = (qJD(2) * t385 - qJD(6)) * t392 + t383 * t409;
	t388 = sin(pkin(10));
	t389 = sin(pkin(6));
	t419 = t388 * t389;
	t390 = cos(pkin(10));
	t418 = t389 * t390;
	t417 = t389 * t392;
	t416 = t389 * t393;
	t391 = cos(pkin(6));
	t415 = t391 * t392;
	t414 = t391 * t393;
	t413 = qJD(2) * t392;
	t412 = qJD(2) * t393;
	t411 = qJD(3) * t383;
	t410 = qJD(3) * t385;
	t386 = pkin(12) + qJ(6);
	t382 = sin(t386);
	t408 = qJD(6) * t382;
	t384 = cos(t386);
	t407 = qJD(6) * t384;
	t406 = qJD(6) * t385;
	t405 = t388 * t415;
	t404 = t389 * t413;
	t403 = t389 * t412;
	t396 = -t388 * t392 + t390 * t414;
	t372 = t396 * qJD(2);
	t400 = -t396 * t406 + t372;
	t378 = t388 * t414 + t390 * t392;
	t374 = t378 * qJD(2);
	t399 = t378 * t406 - t374;
	t398 = (qJD(2) - t406) * t393;
	t377 = t388 * t393 + t390 * t415;
	t366 = -t377 * t383 - t385 * t418;
	t397 = -t377 * t385 + t383 * t418;
	t379 = t390 * t393 - t405;
	t368 = -t379 * t383 + t385 * t419;
	t369 = t379 * t385 + t383 * t419;
	t371 = t391 * t383 + t385 * t417;
	t370 = -t383 * t417 + t391 * t385;
	t373 = t377 * qJD(2);
	t395 = qJD(6) * t377 - t373 * t385 - t396 * t411;
	t375 = -qJD(2) * t405 + t390 * t412;
	t394 = qJD(6) * t379 - t375 * t385 + t378 * t411;
	t365 = t370 * qJD(3) + t385 * t403;
	t364 = -t371 * qJD(3) - t383 * t403;
	t363 = t368 * qJD(3) - t374 * t385;
	t362 = -t369 * qJD(3) + t374 * t383;
	t361 = t366 * qJD(3) + t372 * t385;
	t360 = t397 * qJD(3) - t372 * t383;
	t1 = [0, t399 * t382 + t394 * t384, t362 * t384 - t368 * t408, 0, 0, -t363 * t382 + t375 * t384 + (-t369 * t384 - t378 * t382) * qJD(6); 0, t400 * t382 + t395 * t384, t360 * t384 - t366 * t408, 0, 0, -t361 * t382 + t373 * t384 + (t382 * t396 + t384 * t397) * qJD(6); 0, (t382 * t398 - t384 * t420) * t389, t364 * t384 - t370 * t408, 0, 0, t384 * t404 - t365 * t382 + (-t371 * t384 + t382 * t416) * qJD(6); 0, -t394 * t382 + t399 * t384, -t362 * t382 - t368 * t407, 0, 0, -t363 * t384 - t375 * t382 + (t369 * t382 - t378 * t384) * qJD(6); 0, -t395 * t382 + t400 * t384, -t360 * t382 - t366 * t407, 0, 0, -t361 * t384 - t373 * t382 + (-t382 * t397 + t384 * t396) * qJD(6); 0, (t382 * t420 + t384 * t398) * t389, -t364 * t382 - t370 * t407, 0, 0, -t382 * t404 - t365 * t384 + (t371 * t382 + t384 * t416) * qJD(6); 0, -t375 * t383 - t378 * t410, t363, 0, 0, 0; 0, -t373 * t383 + t396 * t410, t361, 0, 0, 0; 0, (-t383 * t413 + t385 * t409) * t389, t365, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end