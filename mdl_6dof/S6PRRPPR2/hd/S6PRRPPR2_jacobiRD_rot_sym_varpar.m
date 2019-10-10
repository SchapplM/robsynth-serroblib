% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:09
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRRPPR2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPPR2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobiRD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:06
	% EndTime: 2019-10-09 22:09:06
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
	% StartTime: 2019-10-09 22:09:07
	% EndTime: 2019-10-09 22:09:07
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
	% StartTime: 2019-10-09 22:09:07
	% EndTime: 2019-10-09 22:09:07
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
	% StartTime: 2019-10-09 22:09:08
	% EndTime: 2019-10-09 22:09:08
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (66->23), mult. (140->60), div. (0->0), fcn. (148->8), ass. (0->28)
	t279 = sin(pkin(10));
	t280 = sin(pkin(6));
	t295 = t279 * t280;
	t281 = cos(pkin(10));
	t294 = t280 * t281;
	t283 = sin(qJ(2));
	t293 = t280 * t283;
	t282 = cos(pkin(6));
	t292 = t282 * t283;
	t284 = cos(qJ(2));
	t291 = t282 * t284;
	t290 = qJD(2) * t283;
	t278 = qJ(3) + pkin(11);
	t276 = sin(t278);
	t289 = qJD(3) * t276;
	t277 = cos(t278);
	t288 = qJD(3) * t277;
	t287 = qJD(3) * t284;
	t286 = t280 * qJD(2) * t284;
	t272 = -t279 * t283 + t281 * t291;
	t273 = t279 * t284 + t281 * t292;
	t274 = -t279 * t291 - t281 * t283;
	t285 = t279 * t292 - t281 * t284;
	t271 = t285 * qJD(2);
	t270 = t274 * qJD(2);
	t269 = t273 * qJD(2);
	t268 = t272 * qJD(2);
	t1 = [0, t270, 0, 0, 0, 0; 0, t268, 0, 0, 0, 0; 0, t286, 0, 0, 0, 0; 0, -t271 * t277 + t274 * t289, t270 * t276 + (t276 * t295 - t277 * t285) * qJD(3), 0, 0, 0; 0, t269 * t277 + t272 * t289, t268 * t276 + (t273 * t277 - t276 * t294) * qJD(3), 0, 0, 0; 0, (t276 * t287 + t277 * t290) * t280, t276 * t286 + (t276 * t282 + t277 * t293) * qJD(3), 0, 0, 0; 0, t271 * t276 + t274 * t288, t270 * t277 + (t276 * t285 + t277 * t295) * qJD(3), 0, 0, 0; 0, -t269 * t276 + t272 * t288, t268 * t277 + (-t273 * t276 - t277 * t294) * qJD(3), 0, 0, 0; 0, (-t276 * t290 + t277 * t287) * t280, t277 * t286 + (-t276 * t293 + t277 * t282) * qJD(3), 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:09:08
	% EndTime: 2019-10-09 22:09:09
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (246->63), mult. (516->128), div. (0->0), fcn. (564->10), ass. (0->56)
	t384 = sin(pkin(10));
	t385 = sin(pkin(6));
	t417 = t384 * t385;
	t386 = cos(pkin(10));
	t416 = t385 * t386;
	t389 = sin(qJ(2));
	t415 = t385 * t389;
	t387 = cos(pkin(6));
	t414 = t387 * t389;
	t391 = cos(qJ(2));
	t413 = t387 * t391;
	t388 = sin(qJ(6));
	t412 = t388 * t391;
	t390 = cos(qJ(6));
	t411 = t390 * t391;
	t410 = qJD(2) * t389;
	t409 = qJD(2) * t391;
	t383 = qJ(3) + pkin(11);
	t381 = sin(t383);
	t408 = qJD(3) * t381;
	t382 = cos(t383);
	t407 = qJD(3) * t382;
	t406 = qJD(3) * t391;
	t405 = qJD(6) * t381;
	t404 = qJD(6) * t388;
	t403 = qJD(6) * t390;
	t402 = t384 * t414;
	t401 = t385 * t410;
	t400 = t385 * t409;
	t399 = qJD(2) + t405;
	t395 = -t384 * t389 + t386 * t413;
	t371 = t395 * qJD(2);
	t398 = -t395 * t405 - t371;
	t377 = t384 * t413 + t386 * t389;
	t373 = t377 * qJD(2);
	t397 = t377 * t405 + t373;
	t376 = t384 * t391 + t386 * t414;
	t365 = t376 * t381 + t382 * t416;
	t366 = t376 * t382 - t381 * t416;
	t378 = t386 * t391 - t402;
	t396 = -t378 * t381 + t382 * t417;
	t368 = t378 * t382 + t381 * t417;
	t370 = t387 * t381 + t382 * t415;
	t369 = t381 * t415 - t387 * t382;
	t372 = t376 * qJD(2);
	t394 = -qJD(6) * t376 - t372 * t381 + t395 * t407;
	t374 = -qJD(2) * t402 + t386 * t409;
	t393 = -qJD(6) * t378 - t374 * t381 - t377 * t407;
	t392 = t382 * t406 + (-qJD(2) * t381 - qJD(6)) * t389;
	t364 = -t369 * qJD(3) + t382 * t400;
	t363 = t370 * qJD(3) + t381 * t400;
	t362 = t396 * qJD(3) - t373 * t382;
	t361 = t368 * qJD(3) - t373 * t381;
	t360 = -t365 * qJD(3) + t371 * t382;
	t359 = t366 * qJD(3) + t371 * t381;
	t1 = [0, t393 * t388 - t397 * t390, t362 * t388 + t368 * t403, 0, 0, t361 * t390 - t374 * t388 + (-t377 * t390 + t388 * t396) * qJD(6); 0, t394 * t388 - t398 * t390, t360 * t388 + t366 * t403, 0, 0, t359 * t390 - t372 * t388 + (-t365 * t388 + t390 * t395) * qJD(6); 0, (t392 * t388 + t399 * t411) * t385, t364 * t388 + t370 * t403, 0, 0, -t388 * t401 + t363 * t390 + (-t369 * t388 + t385 * t411) * qJD(6); 0, t397 * t388 + t393 * t390, t362 * t390 - t368 * t404, 0, 0, -t361 * t388 - t374 * t390 + (t377 * t388 + t390 * t396) * qJD(6); 0, t398 * t388 + t394 * t390, t360 * t390 - t366 * t404, 0, 0, -t359 * t388 - t372 * t390 + (-t365 * t390 - t388 * t395) * qJD(6); 0, (t392 * t390 - t399 * t412) * t385, t364 * t390 - t370 * t404, 0, 0, -t390 * t401 - t363 * t388 + (-t369 * t390 - t385 * t412) * qJD(6); 0, -t374 * t382 + t377 * t408, -t361, 0, 0, 0; 0, -t372 * t382 - t395 * t408, -t359, 0, 0, 0; 0, (-t381 * t406 - t382 * t410) * t385, -t363, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end