% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% JRD_rot [9x6]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S6PRPRPR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobiRD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_jacobiRD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobiRD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:31
	% EndTime: 2019-10-09 21:35:31
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
	% StartTime: 2019-10-09 21:35:32
	% EndTime: 2019-10-09 21:35:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->7), mult. (42->23), div. (0->0), fcn. (42->8), ass. (0->14)
	t176 = cos(pkin(6));
	t177 = sin(qJ(2));
	t182 = t176 * t177;
	t178 = cos(qJ(2));
	t181 = t176 * t178;
	t180 = qJD(2) * sin(pkin(6));
	t179 = t177 * t180;
	t175 = cos(pkin(10));
	t174 = cos(pkin(11));
	t172 = sin(pkin(10));
	t171 = sin(pkin(11));
	t170 = (t172 * t182 - t175 * t178) * qJD(2);
	t169 = (-t172 * t178 - t175 * t182) * qJD(2);
	t1 = [0, t170 * t174, 0, 0, 0, 0; 0, t169 * t174, 0, 0, 0, 0; 0, -t174 * t179, 0, 0, 0, 0; 0, -t170 * t171, 0, 0, 0, 0; 0, -t169 * t171, 0, 0, 0, 0; 0, t171 * t179, 0, 0, 0, 0; 0, (-t172 * t181 - t175 * t177) * qJD(2), 0, 0, 0, 0; 0, (-t172 * t177 + t175 * t181) * qJD(2), 0, 0, 0, 0; 0, t178 * t180, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:32
	% EndTime: 2019-10-09 21:35:32
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (66->23), mult. (140->60), div. (0->0), fcn. (148->8), ass. (0->28)
	t224 = sin(pkin(10));
	t225 = sin(pkin(6));
	t240 = t224 * t225;
	t226 = cos(pkin(10));
	t239 = t225 * t226;
	t228 = sin(qJ(2));
	t238 = t225 * t228;
	t227 = cos(pkin(6));
	t237 = t227 * t228;
	t229 = cos(qJ(2));
	t236 = t227 * t229;
	t235 = qJD(2) * t228;
	t223 = pkin(11) + qJ(4);
	t221 = sin(t223);
	t234 = qJD(4) * t221;
	t222 = cos(t223);
	t233 = qJD(4) * t222;
	t232 = qJD(4) * t229;
	t231 = t225 * qJD(2) * t229;
	t217 = -t224 * t228 + t226 * t236;
	t218 = t224 * t229 + t226 * t237;
	t219 = -t224 * t236 - t226 * t228;
	t230 = t224 * t237 - t226 * t229;
	t216 = t230 * qJD(2);
	t215 = t219 * qJD(2);
	t214 = t218 * qJD(2);
	t213 = t217 * qJD(2);
	t1 = [0, t216 * t222 - t219 * t234, 0, -t215 * t221 + (-t221 * t240 + t222 * t230) * qJD(4), 0, 0; 0, -t214 * t222 - t217 * t234, 0, -t213 * t221 + (-t218 * t222 + t221 * t239) * qJD(4), 0, 0; 0, (-t221 * t232 - t222 * t235) * t225, 0, -t221 * t231 + (-t221 * t227 - t222 * t238) * qJD(4), 0, 0; 0, -t216 * t221 - t219 * t233, 0, -t215 * t222 + (-t221 * t230 - t222 * t240) * qJD(4), 0, 0; 0, t214 * t221 - t217 * t233, 0, -t213 * t222 + (t218 * t221 + t222 * t239) * qJD(4), 0, 0; 0, (t221 * t235 - t222 * t232) * t225, 0, -t222 * t231 + (t221 * t238 - t222 * t227) * qJD(4), 0, 0; 0, t215, 0, 0, 0, 0; 0, t213, 0, 0, 0, 0; 0, t231, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:33
	% EndTime: 2019-10-09 21:35:33
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (109->33), mult. (240->83), div. (0->0), fcn. (252->10), ass. (0->36)
	t316 = pkin(11) + qJ(4);
	t315 = cos(t316);
	t323 = sin(qJ(2));
	t338 = t315 * t323;
	t318 = sin(pkin(10));
	t319 = sin(pkin(6));
	t337 = t318 * t319;
	t321 = cos(pkin(10));
	t336 = t319 * t321;
	t335 = t319 * t323;
	t322 = cos(pkin(6));
	t334 = t322 * t323;
	t324 = cos(qJ(2));
	t333 = t322 * t324;
	t314 = sin(t316);
	t332 = qJD(4) * t314;
	t331 = qJD(4) * t315;
	t330 = qJD(4) * t324;
	t329 = qJD(2) * t319 * t324;
	t328 = t314 * t330;
	t310 = -t318 * t323 + t321 * t333;
	t311 = t318 * t324 + t321 * t334;
	t312 = -t318 * t333 - t321 * t323;
	t327 = t318 * t334 - t321 * t324;
	t307 = t311 * qJD(2);
	t326 = t307 * t315 + t310 * t332;
	t309 = t327 * qJD(2);
	t325 = -t309 * t315 + t312 * t332;
	t320 = cos(pkin(12));
	t317 = sin(pkin(12));
	t308 = t312 * qJD(2);
	t306 = t310 * qJD(2);
	t305 = -t314 * t329 + (-t314 * t322 - t315 * t335) * qJD(4);
	t304 = -t308 * t314 + (-t314 * t337 + t315 * t327) * qJD(4);
	t303 = -t306 * t314 + (-t311 * t315 + t314 * t336) * qJD(4);
	t1 = [0, t308 * t317 - t325 * t320, 0, t304 * t320, 0, 0; 0, t306 * t317 - t326 * t320, 0, t303 * t320, 0, 0; 0, (-t320 * t328 + (t317 * t324 - t320 * t338) * qJD(2)) * t319, 0, t305 * t320, 0, 0; 0, t308 * t320 + t325 * t317, 0, -t304 * t317, 0, 0; 0, t306 * t320 + t326 * t317, 0, -t303 * t317, 0, 0; 0, (t317 * t328 + (t317 * t338 + t320 * t324) * qJD(2)) * t319, 0, -t305 * t317, 0, 0; 0, t309 * t314 + t312 * t331, 0, t308 * t315 + (t314 * t327 + t315 * t337) * qJD(4), 0, 0; 0, -t307 * t314 + t310 * t331, 0, t306 * t315 + (-t311 * t314 - t315 * t336) * qJD(4), 0, 0; 0, (-qJD(2) * t314 * t323 + t315 * t330) * t319, 0, t315 * t329 + (-t314 * t335 + t315 * t322) * qJD(4), 0, 0;];
	JRD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiRD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:35:33
	% EndTime: 2019-10-09 21:35:34
	% DurationCPUTime: 0.42s
	% Computational Cost: add. (312->61), mult. (516->128), div. (0->0), fcn. (564->10), ass. (0->56)
	t381 = pkin(11) + qJ(4);
	t377 = sin(t381);
	t379 = cos(t381);
	t386 = sin(qJ(2));
	t387 = cos(qJ(2));
	t403 = qJD(4) * t387;
	t414 = (qJD(2) * t379 - qJD(6)) * t386 + t377 * t403;
	t382 = sin(pkin(10));
	t383 = sin(pkin(6));
	t413 = t382 * t383;
	t384 = cos(pkin(10));
	t412 = t383 * t384;
	t411 = t383 * t386;
	t410 = t383 * t387;
	t385 = cos(pkin(6));
	t409 = t385 * t386;
	t408 = t385 * t387;
	t407 = qJD(2) * t386;
	t406 = qJD(2) * t387;
	t405 = qJD(4) * t377;
	t404 = qJD(4) * t379;
	t380 = pkin(12) + qJ(6);
	t376 = sin(t380);
	t402 = qJD(6) * t376;
	t378 = cos(t380);
	t401 = qJD(6) * t378;
	t400 = qJD(6) * t379;
	t399 = t382 * t409;
	t398 = t383 * t407;
	t397 = t383 * t406;
	t390 = -t382 * t386 + t384 * t408;
	t366 = t390 * qJD(2);
	t394 = -t390 * t400 + t366;
	t372 = t382 * t408 + t384 * t386;
	t368 = t372 * qJD(2);
	t393 = t372 * t400 - t368;
	t392 = (qJD(2) - t400) * t387;
	t371 = t382 * t387 + t384 * t409;
	t360 = -t371 * t377 - t379 * t412;
	t391 = -t371 * t379 + t377 * t412;
	t373 = t384 * t387 - t399;
	t362 = -t373 * t377 + t379 * t413;
	t363 = t373 * t379 + t377 * t413;
	t365 = t385 * t377 + t379 * t411;
	t364 = -t377 * t411 + t385 * t379;
	t367 = t371 * qJD(2);
	t389 = qJD(6) * t371 - t367 * t379 - t390 * t405;
	t369 = -qJD(2) * t399 + t384 * t406;
	t388 = qJD(6) * t373 - t369 * t379 + t372 * t405;
	t359 = t364 * qJD(4) + t379 * t397;
	t358 = -t365 * qJD(4) - t377 * t397;
	t357 = t362 * qJD(4) - t368 * t379;
	t356 = -t363 * qJD(4) + t368 * t377;
	t355 = t360 * qJD(4) + t366 * t379;
	t354 = t391 * qJD(4) - t366 * t377;
	t1 = [0, t393 * t376 + t388 * t378, 0, t356 * t378 - t362 * t402, 0, -t357 * t376 + t369 * t378 + (-t363 * t378 - t372 * t376) * qJD(6); 0, t394 * t376 + t389 * t378, 0, t354 * t378 - t360 * t402, 0, -t355 * t376 + t367 * t378 + (t376 * t390 + t378 * t391) * qJD(6); 0, (t376 * t392 - t414 * t378) * t383, 0, t358 * t378 - t364 * t402, 0, t378 * t398 - t359 * t376 + (-t365 * t378 + t376 * t410) * qJD(6); 0, -t388 * t376 + t393 * t378, 0, -t356 * t376 - t362 * t401, 0, -t357 * t378 - t369 * t376 + (t363 * t376 - t372 * t378) * qJD(6); 0, -t389 * t376 + t394 * t378, 0, -t354 * t376 - t360 * t401, 0, -t355 * t378 - t367 * t376 + (-t376 * t391 + t378 * t390) * qJD(6); 0, (t414 * t376 + t378 * t392) * t383, 0, -t358 * t376 - t364 * t401, 0, -t376 * t398 - t359 * t378 + (t365 * t376 + t378 * t410) * qJD(6); 0, -t369 * t377 - t372 * t404, 0, t357, 0, 0; 0, -t367 * t377 + t390 * t404, 0, t355, 0, 0; 0, (-t377 * t407 + t379 * t403) * t383, 0, t359, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,6);
end