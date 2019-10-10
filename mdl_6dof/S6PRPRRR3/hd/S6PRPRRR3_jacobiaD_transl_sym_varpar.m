% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:57
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR3_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR3_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR3_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:34
	% EndTime: 2019-10-09 21:57:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:35
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
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:35
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (29->14), mult. (103->28), div. (0->0), fcn. (90->8), ass. (0->16)
	t153 = sin(pkin(11));
	t156 = cos(pkin(11));
	t159 = cos(qJ(2));
	t157 = cos(pkin(6));
	t158 = sin(qJ(2));
	t165 = t157 * t158;
	t168 = t153 * t165 - t156 * t159;
	t167 = r_i_i_C(3) + qJ(3);
	t164 = t157 * t159;
	t162 = qJD(2) * t167;
	t161 = -cos(pkin(12)) * r_i_i_C(1) + sin(pkin(12)) * r_i_i_C(2) - pkin(2);
	t160 = t153 * t159 + t156 * t165;
	t154 = sin(pkin(6));
	t150 = t168 * qJD(2);
	t148 = t160 * qJD(2);
	t1 = [0, -t168 * qJD(3) - (t153 * t164 + t156 * t158) * t162 - t161 * t150, -t150, 0, 0, 0; 0, t160 * qJD(3) - (t153 * t158 - t156 * t164) * t162 + t161 * t148, t148, 0, 0, 0; 0, (qJD(3) * t158 + (t161 * t158 + t167 * t159) * qJD(2)) * t154, t154 * qJD(2) * t158, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:36
	% DurationCPUTime: 0.27s
	% Computational Cost: add. (98->33), mult. (216->67), div. (0->0), fcn. (199->9), ass. (0->30)
	t212 = r_i_i_C(3) + pkin(8) + qJ(3);
	t192 = sin(pkin(11));
	t193 = sin(pkin(6));
	t211 = t192 * t193;
	t194 = cos(pkin(11));
	t210 = t193 * t194;
	t197 = sin(qJ(2));
	t209 = t193 * t197;
	t195 = cos(pkin(6));
	t208 = t195 * t197;
	t198 = cos(qJ(2));
	t207 = t195 * t198;
	t206 = qJD(2) * t197;
	t205 = qJD(2) * t198;
	t204 = t192 * t206;
	t203 = t194 * t205;
	t191 = pkin(12) + qJ(4);
	t189 = sin(t191);
	t190 = cos(t191);
	t202 = t189 * r_i_i_C(1) + t190 * r_i_i_C(2);
	t201 = -t190 * r_i_i_C(1) + t189 * r_i_i_C(2) - cos(pkin(12)) * pkin(3) - pkin(2);
	t183 = t192 * t198 + t194 * t208;
	t200 = t192 * t207 + t194 * t197;
	t199 = qJD(4) * t202;
	t185 = -t192 * t208 + t194 * t198;
	t181 = -t195 * t204 + t203;
	t180 = t200 * qJD(2);
	t179 = t183 * qJD(2);
	t178 = -t195 * t203 + t204;
	t1 = [0, t185 * qJD(3) - t212 * t180 + t201 * t181 + t200 * t199, t181, t202 * t180 + ((-t185 * t190 - t189 * t211) * r_i_i_C(1) + (t185 * t189 - t190 * t211) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t183 * qJD(3) - t212 * t178 - (-t192 * t197 + t194 * t207) * t199 + t201 * t179, t179, t202 * t178 + ((-t183 * t190 + t189 * t210) * r_i_i_C(1) + (t183 * t189 + t190 * t210) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (qJD(3) * t197 - t198 * t199 + (t201 * t197 + t212 * t198) * qJD(2)) * t193, t193 * t206, -t202 * t193 * t205 + ((-t189 * t195 - t190 * t209) * r_i_i_C(1) + (t189 * t209 - t190 * t195) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:35
	% EndTime: 2019-10-09 21:57:36
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (256->46), mult. (358->85), div. (0->0), fcn. (331->11), ass. (0->43)
	t258 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(3);
	t228 = pkin(12) + qJ(4);
	t226 = qJ(5) + t228;
	t222 = sin(t226);
	t229 = qJD(4) + qJD(5);
	t257 = t222 * t229;
	t223 = cos(t226);
	t256 = t223 * t229;
	t230 = sin(pkin(11));
	t231 = sin(pkin(6));
	t255 = t230 * t231;
	t232 = cos(pkin(11));
	t254 = t231 * t232;
	t234 = sin(qJ(2));
	t253 = t231 * t234;
	t233 = cos(pkin(6));
	t252 = t233 * t234;
	t235 = cos(qJ(2));
	t251 = t233 * t235;
	t216 = t230 * t235 + t232 * t252;
	t246 = qJD(2) * t235;
	t242 = t232 * t246;
	t247 = qJD(2) * t234;
	t243 = t230 * t247;
	t211 = -t233 * t242 + t243;
	t241 = t229 * t254 + t211;
	t250 = (-t216 * t256 + t241 * t222) * r_i_i_C(1) + (t216 * t257 + t241 * t223) * r_i_i_C(2);
	t218 = -t230 * t252 + t232 * t235;
	t238 = t230 * t251 + t232 * t234;
	t213 = t238 * qJD(2);
	t240 = -t229 * t255 + t213;
	t249 = (-t218 * t256 + t240 * t222) * r_i_i_C(1) + (t218 * t257 + t240 * t223) * r_i_i_C(2);
	t244 = t231 * t246;
	t237 = -t229 * t233 - t244;
	t245 = t229 * t253;
	t248 = (t237 * t222 - t223 * t245) * r_i_i_C(1) + (t222 * t245 + t237 * t223) * r_i_i_C(2);
	t225 = cos(t228);
	t239 = -r_i_i_C(1) * t223 + r_i_i_C(2) * t222 - pkin(4) * t225 - cos(pkin(12)) * pkin(3) - pkin(2);
	t224 = sin(t228);
	t236 = -pkin(4) * qJD(4) * t224 + (-r_i_i_C(1) * t222 - r_i_i_C(2) * t223) * t229;
	t214 = -t233 * t243 + t242;
	t212 = t216 * qJD(2);
	t1 = [0, t218 * qJD(3) - t258 * t213 + t239 * t214 - t236 * t238, t214, (t213 * t224 + (-t218 * t225 - t224 * t255) * qJD(4)) * pkin(4) + t249, t249, 0; 0, t216 * qJD(3) - t258 * t211 + t236 * (-t230 * t234 + t232 * t251) + t239 * t212, t212, (t211 * t224 + (-t216 * t225 + t224 * t254) * qJD(4)) * pkin(4) + t250, t250, 0; 0, (qJD(3) * t234 + t236 * t235 + (t239 * t234 + t258 * t235) * qJD(2)) * t231, t231 * t247, (-t224 * t244 + (-t224 * t233 - t225 * t253) * qJD(4)) * pkin(4) + t248, t248, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:57:37
	% EndTime: 2019-10-09 21:57:37
	% DurationCPUTime: 0.70s
	% Computational Cost: add. (801->93), mult. (1063->161), div. (0->0), fcn. (1055->13), ass. (0->62)
	t408 = pkin(10) + r_i_i_C(3);
	t367 = sin(qJ(6));
	t369 = cos(qJ(6));
	t381 = t369 * r_i_i_C(1) - t367 * r_i_i_C(2);
	t414 = pkin(5) + t381;
	t394 = qJD(6) * t369;
	t395 = qJD(6) * t367;
	t413 = -r_i_i_C(1) * t395 - t394 * r_i_i_C(2);
	t363 = pkin(12) + qJ(4);
	t361 = qJ(5) + t363;
	t357 = sin(t361);
	t358 = cos(t361);
	t359 = sin(t363);
	t364 = qJD(4) + qJD(5);
	t412 = -pkin(4) * qJD(4) * t359 - (t357 * t414 - t408 * t358) * t364;
	t368 = sin(qJ(2));
	t370 = cos(qJ(2));
	t365 = sin(pkin(11));
	t404 = cos(pkin(6));
	t387 = t365 * t404;
	t403 = cos(pkin(11));
	t350 = -t368 * t387 + t403 * t370;
	t410 = t367 * r_i_i_C(1) + t369 * r_i_i_C(2);
	t402 = t357 * t364;
	t401 = t358 * t364;
	t366 = sin(pkin(6));
	t400 = t365 * t366;
	t399 = t366 * t368;
	t398 = t366 * t370;
	t397 = qJD(2) * t368;
	t396 = qJD(6) * t358;
	t391 = t357 * t399;
	t390 = t358 * t399;
	t389 = qJD(2) * t398;
	t388 = t366 * t397;
	t386 = t366 * t403;
	t382 = t358 * t386;
	t349 = t403 * t368 + t370 * t387;
	t345 = t349 * qJD(2);
	t380 = t364 * t400 - t345;
	t379 = t404 * t403;
	t377 = t370 * t379;
	t376 = t404 * t364 + t389;
	t348 = t365 * t370 + t368 * t379;
	t360 = cos(t363);
	t375 = -t408 * t357 - t414 * t358 - pkin(4) * t360 - cos(pkin(12)) * pkin(3) - pkin(2);
	t343 = -qJD(2) * t377 + t365 * t397;
	t327 = -t343 * t358 - t348 * t402 - t364 * t382;
	t374 = t413 * (-t348 * t357 - t382) + t408 * t327 + t414 * (-t348 * t401 + (t364 * t386 + t343) * t357);
	t329 = -t350 * t402 + t380 * t358;
	t373 = t413 * (-t350 * t357 + t358 * t400) + t408 * t329 + t414 * (-t350 * t401 - t380 * t357);
	t334 = t376 * t358 - t364 * t391;
	t372 = t413 * (t404 * t358 - t391) + t408 * t334 + t414 * (-t376 * t357 - t364 * t390);
	t371 = t410 * t396 - t412;
	t362 = -pkin(9) - pkin(8) - qJ(3);
	t347 = t365 * t368 - t377;
	t346 = t350 * qJD(2);
	t344 = t348 * qJD(2);
	t342 = t404 * t357 + t390;
	t338 = t350 * t358 + t357 * t400;
	t336 = t348 * t358 - t357 * t386;
	t1 = [0, (-t345 * t367 + t350 * t394) * r_i_i_C(1) + (-t345 * t369 - t350 * t395) * r_i_i_C(2) + t345 * t362 + t350 * qJD(3) + t375 * t346 + t371 * t349, t346, (t345 * t359 + (-t350 * t360 - t359 * t400) * qJD(4)) * pkin(4) + t373, t373, (-t329 * t367 + t346 * t369) * r_i_i_C(1) + (-t329 * t369 - t346 * t367) * r_i_i_C(2) + ((-t338 * t369 - t349 * t367) * r_i_i_C(1) + (t338 * t367 - t349 * t369) * r_i_i_C(2)) * qJD(6); 0, (-t343 * t367 + t348 * t394) * r_i_i_C(1) + (-t343 * t369 - t348 * t395) * r_i_i_C(2) + t343 * t362 + t348 * qJD(3) + t375 * t344 + t371 * t347, t344, (t343 * t359 + (-t348 * t360 + t359 * t386) * qJD(4)) * pkin(4) + t374, t374, (-t327 * t367 + t344 * t369) * r_i_i_C(1) + (-t327 * t369 - t344 * t367) * r_i_i_C(2) + ((-t336 * t369 - t347 * t367) * r_i_i_C(1) + (t336 * t367 - t347 * t369) * r_i_i_C(2)) * qJD(6); 0, ((t375 * qJD(2) + t381 * qJD(6) + qJD(3)) * t368 + (-qJD(2) * t362 + t410 * (qJD(2) - t396) + t412) * t370) * t366, t388, (-t359 * t389 + (-t404 * t359 - t360 * t399) * qJD(4)) * pkin(4) + t372, t372, (-t334 * t367 + t369 * t388) * r_i_i_C(1) + (-t334 * t369 - t367 * t388) * r_i_i_C(2) + ((-t342 * t369 + t367 * t398) * r_i_i_C(1) + (t342 * t367 + t369 * t398) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end