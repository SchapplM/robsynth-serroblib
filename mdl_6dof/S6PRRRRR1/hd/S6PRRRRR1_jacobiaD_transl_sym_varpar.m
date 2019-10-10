% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(12));
	t50 = sin(pkin(12));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (55->26), mult. (191->57), div. (0->0), fcn. (172->8), ass. (0->24)
	t180 = sin(pkin(12));
	t182 = cos(pkin(12));
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
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (180->39), mult. (333->75), div. (0->0), fcn. (304->10), ass. (0->38)
	t224 = sin(pkin(12));
	t226 = cos(pkin(12));
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
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:46
	% DurationCPUTime: 0.22s
	% Computational Cost: add. (411->54), mult. (480->90), div. (0->0), fcn. (433->12), ass. (0->49)
	t236 = sin(pkin(12));
	t238 = cos(pkin(12));
	t241 = sin(qJ(2));
	t239 = cos(pkin(6));
	t243 = cos(qJ(2));
	t259 = t239 * t243;
	t271 = -t236 * t241 + t238 * t259;
	t235 = qJ(3) + qJ(4);
	t230 = sin(t235);
	t270 = pkin(4) * t230;
	t269 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
	t268 = pkin(3) * qJD(3);
	t232 = qJ(5) + t235;
	t227 = sin(t232);
	t233 = qJD(3) + qJD(4);
	t229 = qJD(5) + t233;
	t267 = t227 * t229;
	t228 = cos(t232);
	t266 = t228 * t229;
	t231 = cos(t235);
	t265 = t231 * t233;
	t237 = sin(pkin(6));
	t264 = t236 * t237;
	t262 = t237 * t238;
	t261 = t237 * t241;
	t260 = t239 * t241;
	t219 = t236 * t243 + t238 * t260;
	t214 = t271 * qJD(2);
	t251 = t229 * t262 - t214;
	t258 = (-t219 * t266 + t251 * t227) * r_i_i_C(1) + (t219 * t267 + t251 * t228) * r_i_i_C(2);
	t247 = t236 * t260 - t238 * t243;
	t248 = t236 * t259 + t238 * t241;
	t216 = t248 * qJD(2);
	t250 = -t229 * t264 + t216;
	t257 = (t250 * t227 + t247 * t266) * r_i_i_C(1) + (t250 * t228 - t247 * t267) * r_i_i_C(2);
	t255 = qJD(2) * t243;
	t252 = t237 * t255;
	t246 = -t229 * t239 - t252;
	t254 = t229 * t261;
	t256 = (t246 * t227 - t228 * t254) * r_i_i_C(1) + (t227 * t254 + t246 * t228) * r_i_i_C(2);
	t242 = cos(qJ(3));
	t249 = t242 * pkin(3) + pkin(4) * t231 + r_i_i_C(1) * t228 - r_i_i_C(2) * t227 + pkin(2);
	t245 = qJD(2) * t249;
	t240 = sin(qJ(3));
	t222 = -t233 * t270 - t240 * t268;
	t244 = t222 + (-r_i_i_C(1) * t227 - r_i_i_C(2) * t228) * t229;
	t225 = -t240 * pkin(3) - t270;
	t223 = -pkin(4) * t265 - t242 * t268;
	t1 = [0, -t269 * t216 - t244 * t248 + t247 * t245, -t216 * t225 + t222 * t264 - t223 * t247 + t257, (t247 * t265 + (-t233 * t264 + t216) * t230) * pkin(4) + t257, t257, 0; 0, t269 * t214 - t219 * t245 + t244 * t271, t214 * t225 + t219 * t223 - t222 * t262 + t258, (-t219 * t265 + (t233 * t262 - t214) * t230) * pkin(4) + t258, t258, 0; 0, (t244 * t243 + (-t249 * t241 + t269 * t243) * qJD(2)) * t237, t239 * t222 + (t223 * t241 + t225 * t255) * t237 + t256, (-t261 * t265 + (-t233 * t239 - t252) * t230) * pkin(4) + t256, t256, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:47
	% EndTime: 2019-10-09 23:13:48
	% DurationCPUTime: 0.75s
	% Computational Cost: add. (1182->103), mult. (1356->173), div. (0->0), fcn. (1329->14), ass. (0->72)
	t423 = pkin(11) + r_i_i_C(3);
	t377 = sin(qJ(6));
	t380 = cos(qJ(6));
	t393 = r_i_i_C(1) * t380 - r_i_i_C(2) * t377;
	t429 = pkin(5) + t393;
	t405 = qJD(6) * t380;
	t406 = qJD(6) * t377;
	t428 = -r_i_i_C(1) * t406 - t405 * r_i_i_C(2);
	t372 = qJD(3) + qJD(4);
	t378 = sin(qJ(3));
	t418 = pkin(3) * qJD(3);
	t374 = qJ(3) + qJ(4);
	t369 = sin(t374);
	t422 = pkin(4) * t369;
	t359 = -t372 * t422 - t378 * t418;
	t371 = qJ(5) + t374;
	t366 = sin(t371);
	t367 = cos(t371);
	t368 = qJD(5) + t372;
	t427 = -(t366 * t429 - t423 * t367) * t368 + t359;
	t379 = sin(qJ(2));
	t382 = cos(qJ(2));
	t375 = sin(pkin(12));
	t417 = cos(pkin(6));
	t399 = t375 * t417;
	t416 = cos(pkin(12));
	t358 = -t379 * t399 + t416 * t382;
	t425 = r_i_i_C(1) * t377 + r_i_i_C(2) * t380;
	t415 = t366 * t368;
	t414 = t367 * t368;
	t370 = cos(t374);
	t413 = t370 * t372;
	t376 = sin(pkin(6));
	t412 = t375 * t376;
	t411 = t376 * t379;
	t410 = t376 * t382;
	t409 = qJD(2) * t379;
	t408 = qJD(2) * t382;
	t407 = qJD(6) * t367;
	t403 = t366 * t411;
	t402 = t367 * t411;
	t401 = t376 * t409;
	t400 = t376 * t408;
	t398 = t376 * t416;
	t394 = t367 * t398;
	t357 = t416 * t379 + t382 * t399;
	t351 = t357 * qJD(2);
	t392 = t368 * t412 - t351;
	t391 = t417 * t416;
	t389 = t382 * t391;
	t388 = t417 * t368 + t400;
	t356 = t375 * t382 + t379 * t391;
	t381 = cos(qJ(3));
	t387 = -pkin(3) * t381 - pkin(4) * t370 - t423 * t366 - t367 * t429 - pkin(2);
	t335 = -t358 * t415 + t392 * t367;
	t386 = t428 * (-t358 * t366 + t367 * t412) + t423 * t335 + t429 * (-t358 * t414 - t392 * t366);
	t349 = -qJD(2) * t389 + t375 * t409;
	t333 = -t349 * t367 - t356 * t415 - t368 * t394;
	t385 = t428 * (-t356 * t366 - t394) + t423 * t333 + t429 * (-t356 * t414 + (t368 * t398 + t349) * t366);
	t340 = t388 * t367 - t368 * t403;
	t384 = t428 * (t417 * t367 - t403) + t423 * t340 + t429 * (-t388 * t366 - t368 * t402);
	t383 = t425 * t407 - t427;
	t373 = -pkin(10) - pkin(9) - pkin(8);
	t362 = -pkin(3) * t378 - t422;
	t360 = -pkin(4) * t413 - t381 * t418;
	t355 = t375 * t379 - t389;
	t352 = t358 * qJD(2);
	t350 = t356 * qJD(2);
	t348 = t417 * t366 + t402;
	t344 = t358 * t367 + t366 * t412;
	t342 = t356 * t367 - t366 * t398;
	t1 = [0, (-t351 * t377 + t358 * t405) * r_i_i_C(1) + (-t351 * t380 - t358 * t406) * r_i_i_C(2) + t351 * t373 + t387 * t352 + t383 * t357, -t351 * t362 + t358 * t360 + t359 * t412 + t386, (-t358 * t413 + (-t372 * t412 + t351) * t369) * pkin(4) + t386, t386, (-t335 * t377 + t352 * t380) * r_i_i_C(1) + (-t335 * t380 - t352 * t377) * r_i_i_C(2) + ((-t344 * t380 - t357 * t377) * r_i_i_C(1) + (t344 * t377 - t357 * t380) * r_i_i_C(2)) * qJD(6); 0, (-t349 * t377 + t356 * t405) * r_i_i_C(1) + (-t349 * t380 - t356 * t406) * r_i_i_C(2) + t349 * t373 + t387 * t350 + t383 * t355, -t349 * t362 + t356 * t360 - t359 * t398 + t385, (-t356 * t413 + (t372 * t398 + t349) * t369) * pkin(4) + t385, t385, (-t333 * t377 + t350 * t380) * r_i_i_C(1) + (-t333 * t380 - t350 * t377) * r_i_i_C(2) + ((-t342 * t380 - t355 * t377) * r_i_i_C(1) + (t342 * t377 - t355 * t380) * r_i_i_C(2)) * qJD(6); 0, ((t387 * qJD(2) + t393 * qJD(6)) * t379 + (-qJD(2) * t373 + t425 * (qJD(2) - t407) + t427) * t382) * t376, t417 * t359 + (t360 * t379 + t362 * t408) * t376 + t384, (-t411 * t413 + (-t417 * t372 - t400) * t369) * pkin(4) + t384, t384, (-t340 * t377 + t380 * t401) * r_i_i_C(1) + (-t340 * t380 - t377 * t401) * r_i_i_C(2) + ((-t348 * t380 + t377 * t410) * r_i_i_C(1) + (t348 * t377 + t380 * t410) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end