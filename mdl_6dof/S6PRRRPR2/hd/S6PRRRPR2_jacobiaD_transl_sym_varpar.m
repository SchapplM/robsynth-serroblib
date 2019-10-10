% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRPR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(11));
	t50 = sin(pkin(11));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(11));
	t182 = cos(pkin(11));
	t185 = sin(qJ(2));
	t183 = cos(pkin(6));
	t187 = cos(qJ(2));
	t195 = t183 * t187;
	t201 = -t180 * t185 + t182 * t195;
	t200 = pkin(8) + r_i_i_C(3);
	t181 = sin(pkin(6));
	t184 = sin(qJ(3));
	t198 = t181 * t184;
	t186 = cos(qJ(3));
	t197 = t181 * t186;
	t196 = t183 * t185;
	t193 = t184 * r_i_i_C(1) + t186 * r_i_i_C(2);
	t192 = t186 * r_i_i_C(1) - t184 * r_i_i_C(2) + pkin(2);
	t176 = t180 * t187 + t182 * t196;
	t191 = t180 * t195 + t182 * t185;
	t190 = t180 * t196 - t182 * t187;
	t189 = qJD(3) * t193;
	t188 = qJD(2) * t192;
	t173 = t191 * qJD(2);
	t171 = t201 * qJD(2);
	t1 = [0, -t200 * t173 + t190 * t188 + t191 * t189, t193 * t173 + ((-t180 * t198 + t186 * t190) * r_i_i_C(1) + (-t180 * t197 - t184 * t190) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, t200 * t171 - t176 * t188 - t201 * t189, -t193 * t171 + ((-t176 * t186 + t182 * t198) * r_i_i_C(1) + (t176 * t184 + t182 * t197) * r_i_i_C(2)) * qJD(3), 0, 0, 0; 0, (-t187 * t189 + (-t192 * t185 + t200 * t187) * qJD(2)) * t181, -t193 * t187 * t181 * qJD(2) + ((-t183 * t184 - t185 * t197) * r_i_i_C(1) + (-t183 * t186 + t185 * t198) * r_i_i_C(2)) * qJD(3), 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (180->39), mult. (333->75), div. (0->0), fcn. (304->10), ass. (0->38)
	t224 = sin(pkin(11));
	t226 = cos(pkin(11));
	t229 = sin(qJ(2));
	t227 = cos(pkin(6));
	t231 = cos(qJ(2));
	t247 = t227 * t231;
	t256 = -t224 * t229 + t226 * t247;
	t255 = r_i_i_C(3) + pkin(9) + pkin(8);
	t223 = qJ(3) + qJ(4);
	t220 = sin(t223);
	t222 = qJD(3) + qJD(4);
	t254 = t220 * t222;
	t221 = cos(t223);
	t253 = t221 * t222;
	t225 = sin(pkin(6));
	t252 = t222 * t225;
	t228 = sin(qJ(3));
	t250 = t225 * t228;
	t249 = t225 * t229;
	t248 = t227 * t229;
	t215 = t224 * t231 + t226 * t248;
	t210 = t256 * qJD(2);
	t240 = t226 * t252 - t210;
	t246 = (-t215 * t253 + t240 * t220) * r_i_i_C(1) + (t215 * t254 + t240 * t221) * r_i_i_C(2);
	t236 = t224 * t248 - t226 * t231;
	t237 = t224 * t247 + t226 * t229;
	t212 = t237 * qJD(2);
	t239 = -t224 * t252 + t212;
	t245 = (t239 * t220 + t236 * t253) * r_i_i_C(1) + (t239 * t221 - t236 * t254) * r_i_i_C(2);
	t241 = qJD(2) * t225 * t231;
	t235 = -t222 * t227 - t241;
	t243 = t222 * t249;
	t244 = (t235 * t220 - t221 * t243) * r_i_i_C(1) + (t220 * t243 + t235 * t221) * r_i_i_C(2);
	t230 = cos(qJ(3));
	t238 = t230 * pkin(3) + r_i_i_C(1) * t221 - r_i_i_C(2) * t220 + pkin(2);
	t234 = qJD(2) * t238;
	t233 = -pkin(3) * qJD(3) * t228 + (-r_i_i_C(1) * t220 - r_i_i_C(2) * t221) * t222;
	t1 = [0, -t255 * t212 - t233 * t237 + t236 * t234, (t212 * t228 + (-t224 * t250 + t230 * t236) * qJD(3)) * pkin(3) + t245, t245, 0, 0; 0, t255 * t210 - t215 * t234 + t233 * t256, (-t210 * t228 + (-t215 * t230 + t226 * t250) * qJD(3)) * pkin(3) + t246, t246, 0, 0; 0, (t233 * t231 + (-t238 * t229 + t255 * t231) * qJD(2)) * t225, (-t228 * t241 + (-t227 * t228 - t230 * t249) * qJD(3)) * pkin(3) + t244, t244, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:36
	% EndTime: 2019-10-09 22:48:36
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (471->53), mult. (787->91), div. (0->0), fcn. (763->12), ass. (0->47)
	t362 = r_i_i_C(3) + qJ(5);
	t327 = sin(pkin(12));
	t330 = cos(pkin(12));
	t347 = r_i_i_C(1) * t330 - r_i_i_C(2) * t327 + pkin(4);
	t326 = qJ(3) + qJ(4);
	t323 = sin(t326);
	t324 = cos(t326);
	t335 = cos(qJ(3));
	t365 = t335 * pkin(3) + t362 * t323 + t347 * t324 + pkin(2);
	t325 = qJD(3) + qJD(4);
	t361 = t323 * t325;
	t360 = t324 * t325;
	t328 = sin(pkin(11));
	t329 = sin(pkin(6));
	t359 = t328 * t329;
	t331 = cos(pkin(11));
	t358 = t329 * t331;
	t333 = sin(qJ(3));
	t357 = t329 * t333;
	t334 = sin(qJ(2));
	t356 = t329 * t334;
	t332 = cos(pkin(6));
	t355 = t332 * t334;
	t336 = cos(qJ(2));
	t354 = t332 * t336;
	t353 = qJD(2) * t334;
	t352 = t324 * t356;
	t351 = t323 * t358;
	t350 = t331 * t354;
	t349 = qJD(2) * t329 * t336;
	t345 = t328 * t354 + t331 * t334;
	t312 = t345 * qJD(2);
	t348 = t325 * t359 - t312;
	t346 = t327 * r_i_i_C(1) + t330 * r_i_i_C(2) + pkin(8) + pkin(9);
	t315 = t328 * t336 + t331 * t355;
	t344 = t328 * t355 - t331 * t336;
	t343 = t325 * t332 + t349;
	t310 = -qJD(2) * t350 + t328 * t353;
	t296 = -t310 * t323 + t315 * t360 - t325 * t351;
	t342 = -(-t315 * t324 + t351) * qJD(5) + t362 * (-t310 * t324 - t315 * t361 - t358 * t360) - t347 * t296;
	t298 = t348 * t323 - t344 * t360;
	t341 = -(-t323 * t359 + t324 * t344) * qJD(5) + t362 * (t348 * t324 + t344 * t361) - t347 * t298;
	t303 = t343 * t323 + t325 * t352;
	t340 = -(-t332 * t323 - t352) * qJD(5) + t362 * (t343 * t324 - t356 * t361) - t347 * t303;
	t339 = qJD(2) * t365;
	t338 = -qJD(3) * t333 * pkin(3) + t323 * qJD(5) + (-t347 * t323 + t362 * t324) * t325;
	t1 = [0, -t346 * t312 - t338 * t345 + t344 * t339, (t312 * t333 + (-t328 * t357 + t335 * t344) * qJD(3)) * pkin(3) + t341, t341, t298, 0; 0, -t346 * t310 - t315 * t339 + t338 * (-t328 * t334 + t350), (t310 * t333 + (-t315 * t335 + t331 * t357) * qJD(3)) * pkin(3) + t342, t342, t296, 0; 0, (-t365 * t353 + (t346 * qJD(2) + t338) * t336) * t329, (-t333 * t349 + (-t332 * t333 - t335 * t356) * qJD(3)) * pkin(3) + t340, t340, t303, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:36
	% EndTime: 2019-10-09 22:48:37
	% DurationCPUTime: 0.75s
	% Computational Cost: add. (761->87), mult. (1153->149), div. (0->0), fcn. (1154->14), ass. (0->62)
	t370 = pkin(12) + qJ(6);
	t366 = sin(t370);
	t367 = cos(t370);
	t395 = t366 * r_i_i_C(1) + t367 * r_i_i_C(2);
	t423 = t395 * qJD(6);
	t417 = r_i_i_C(3) + pkin(10) + qJ(5);
	t396 = t367 * r_i_i_C(1) - t366 * r_i_i_C(2);
	t392 = cos(pkin(12)) * pkin(5) + pkin(4) + t396;
	t372 = qJ(3) + qJ(4);
	t368 = sin(t372);
	t369 = cos(t372);
	t371 = qJD(3) + qJD(4);
	t377 = sin(qJ(3));
	t382 = (t392 * t368 - t417 * t369) * t371 + qJD(3) * t377 * pkin(3) + t369 * t423 - t368 * qJD(5);
	t378 = sin(qJ(2));
	t380 = cos(qJ(2));
	t374 = sin(pkin(11));
	t416 = cos(pkin(6));
	t401 = t374 * t416;
	t415 = cos(pkin(11));
	t358 = -t378 * t401 + t415 * t380;
	t414 = t368 * t371;
	t413 = t369 * t371;
	t375 = sin(pkin(6));
	t412 = t374 * t375;
	t411 = t375 * t378;
	t410 = t375 * t380;
	t409 = qJD(2) * t378;
	t405 = t368 * t411;
	t404 = t369 * t411;
	t403 = t375 * t409;
	t402 = qJD(2) * t410;
	t400 = t375 * t415;
	t397 = t368 * t400;
	t357 = t415 * t378 + t380 * t401;
	t353 = t357 * qJD(2);
	t394 = t371 * t412 - t353;
	t393 = t416 * t415;
	t391 = t380 * t393;
	t389 = t396 * qJD(6);
	t388 = t416 * t371 + t402;
	t387 = sin(pkin(12)) * pkin(5) + pkin(9) + pkin(8) + t395;
	t356 = t374 * t380 + t378 * t393;
	t379 = cos(qJ(3));
	t386 = -t379 * pkin(3) - t417 * t368 - t392 * t369 - pkin(2);
	t351 = -qJD(2) * t391 + t374 * t409;
	t335 = -t351 * t368 + t356 * t413 - t371 * t397;
	t336 = -t356 * t414 + (-t371 * t400 - t351) * t369;
	t346 = t356 * t369 - t397;
	t385 = t346 * qJD(5) - t423 * (-t356 * t368 - t369 * t400) + t417 * t336 - t392 * t335;
	t337 = t358 * t413 + t394 * t368;
	t338 = -t358 * t414 + t394 * t369;
	t348 = t358 * t369 + t368 * t412;
	t384 = t348 * qJD(5) - t423 * (-t358 * t368 + t369 * t412) + t417 * t338 - t392 * t337;
	t343 = t388 * t368 + t371 * t404;
	t344 = t388 * t369 - t371 * t405;
	t350 = t416 * t368 + t404;
	t383 = t350 * qJD(5) - t423 * (t416 * t369 - t405) + t417 * t344 - t392 * t343;
	t355 = t374 * t378 - t391;
	t354 = t358 * qJD(2);
	t352 = t356 * qJD(2);
	t1 = [0, -t387 * t353 + t386 * t354 + t382 * t357 + t358 * t389, (t353 * t377 + (-t358 * t379 - t377 * t412) * qJD(3)) * pkin(3) + t384, t384, t337, (-t338 * t366 + t354 * t367) * r_i_i_C(1) + (-t338 * t367 - t354 * t366) * r_i_i_C(2) + ((-t348 * t367 - t357 * t366) * r_i_i_C(1) + (t348 * t366 - t357 * t367) * r_i_i_C(2)) * qJD(6); 0, -t387 * t351 + t386 * t352 + t382 * t355 + t356 * t389, (t351 * t377 + (-t356 * t379 + t377 * t400) * qJD(3)) * pkin(3) + t385, t385, t335, (-t336 * t366 + t352 * t367) * r_i_i_C(1) + (-t336 * t367 - t352 * t366) * r_i_i_C(2) + ((-t346 * t367 - t355 * t366) * r_i_i_C(1) + (t346 * t366 - t355 * t367) * r_i_i_C(2)) * qJD(6); 0, ((t386 * qJD(2) + t389) * t378 + (t387 * qJD(2) - t382) * t380) * t375, (-t377 * t402 + (-t416 * t377 - t379 * t411) * qJD(3)) * pkin(3) + t383, t383, t343, (-t344 * t366 + t367 * t403) * r_i_i_C(1) + (-t344 * t367 - t366 * t403) * r_i_i_C(2) + ((-t350 * t367 + t366 * t410) * r_i_i_C(1) + (t350 * t366 + t367 * t410) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end