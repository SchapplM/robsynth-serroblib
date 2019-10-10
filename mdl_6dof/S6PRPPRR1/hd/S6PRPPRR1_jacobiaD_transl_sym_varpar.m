% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:24
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPPRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiaD_transl_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->8)
	t53 = cos(pkin(6));
	t54 = sin(qJ(2));
	t57 = t53 * t54;
	t55 = cos(qJ(2));
	t56 = t53 * t55;
	t52 = cos(pkin(10));
	t50 = sin(pkin(10));
	t1 = [0, ((t50 * t57 - t52 * t55) * r_i_i_C(1) + (t50 * t56 + t52 * t54) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, ((-t50 * t55 - t52 * t57) * r_i_i_C(1) + (t50 * t54 - t52 * t56) * r_i_i_C(2)) * qJD(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t54 - r_i_i_C(2) * t55) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (22->14), mult. (83->33), div. (0->0), fcn. (72->8), ass. (0->16)
	t89 = sin(pkin(11));
	t92 = cos(pkin(11));
	t95 = sin(qJ(2));
	t96 = cos(qJ(2));
	t98 = t89 * t96 + t92 * t95;
	t88 = t98 * qJD(2);
	t94 = cos(pkin(6));
	t100 = t94 * t95;
	t99 = pkin(2) * qJD(2);
	t97 = t89 * t95 - t92 * t96;
	t87 = t97 * qJD(2);
	t93 = cos(pkin(10));
	t90 = sin(pkin(10));
	t86 = t94 * t88;
	t85 = t94 * t87;
	t1 = [0, (t90 * t86 + t93 * t87) * r_i_i_C(1) + (-t90 * t85 + t93 * t88) * r_i_i_C(2) + (t90 * t100 - t93 * t96) * t99, 0, 0, 0, 0; 0, (-t93 * t86 + t90 * t87) * r_i_i_C(1) + (t93 * t85 + t90 * t88) * r_i_i_C(2) + (-t93 * t100 - t90 * t96) * t99, 0, 0, 0, 0; 0, (-t95 * pkin(2) - t98 * r_i_i_C(1) + t97 * r_i_i_C(2)) * sin(pkin(6)) * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (69->24), mult. (241->43), div. (0->0), fcn. (228->10), ass. (0->23)
	t203 = -r_i_i_C(3) - qJ(4);
	t202 = qJD(2) * pkin(2);
	t196 = cos(pkin(6));
	t197 = sin(qJ(2));
	t201 = t196 * t197;
	t190 = sin(pkin(11));
	t194 = cos(pkin(11));
	t198 = cos(qJ(2));
	t200 = t190 * t198 + t194 * t197;
	t188 = t197 * t190 - t198 * t194;
	t199 = cos(pkin(12)) * r_i_i_C(1) - sin(pkin(12)) * r_i_i_C(2) + pkin(3);
	t185 = t200 * t196;
	t187 = t200 * qJD(2);
	t186 = t188 * qJD(2);
	t195 = cos(pkin(10));
	t192 = sin(pkin(6));
	t191 = sin(pkin(10));
	t184 = qJD(2) * t185;
	t183 = t196 * t186;
	t182 = t192 * t187;
	t179 = t191 * t184 + t195 * t186;
	t177 = -t195 * t184 + t191 * t186;
	t1 = [0, -(t191 * t185 + t195 * t188) * qJD(4) + t203 * (-t191 * t183 + t195 * t187) + (t191 * t201 - t195 * t198) * t202 + t199 * t179, 0, -t179, 0, 0; 0, -(-t195 * t185 + t191 * t188) * qJD(4) + t203 * (t195 * t183 + t191 * t187) + (-t191 * t198 - t195 * t201) * t202 + t199 * t177, 0, -t177, 0, 0; 0, -t199 * t182 + (t200 * qJD(4) + t203 * t186 - t197 * t202) * t192, 0, t182, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:29
	% EndTime: 2019-10-09 21:24:29
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (168->46), mult. (442->88), div. (0->0), fcn. (445->11), ass. (0->40)
	t259 = cos(pkin(6));
	t254 = sin(pkin(11));
	t257 = cos(pkin(11));
	t261 = sin(qJ(2));
	t262 = cos(qJ(2));
	t265 = t262 * t254 + t261 * t257;
	t240 = t265 * t259;
	t272 = qJD(2) * t262;
	t273 = qJD(2) * t261;
	t279 = t254 * t273 - t257 * t272;
	t278 = -r_i_i_C(3) - pkin(8) - qJ(4);
	t277 = pkin(2) * qJD(2);
	t255 = sin(pkin(10));
	t256 = sin(pkin(6));
	t276 = t255 * t256;
	t258 = cos(pkin(10));
	t275 = t256 * t258;
	t274 = t259 * t261;
	t253 = pkin(12) + qJ(5);
	t251 = sin(t253);
	t252 = cos(t253);
	t269 = t251 * r_i_i_C(1) + t252 * r_i_i_C(2);
	t235 = t279 * t259;
	t242 = -t254 * t272 - t257 * t273;
	t268 = t258 * t235 - t255 * t242;
	t267 = t255 * t235 + t258 * t242;
	t243 = t261 * t254 - t262 * t257;
	t230 = t258 * t240 - t255 * t243;
	t266 = t255 * t240 + t258 * t243;
	t264 = t252 * r_i_i_C(1) - t251 * r_i_i_C(2) + cos(pkin(12)) * pkin(4) + pkin(3);
	t263 = qJD(5) * t269;
	t238 = t265 * t256;
	t241 = t243 * qJD(2);
	t239 = t243 * t259;
	t236 = qJD(2) * t240;
	t234 = qJD(2) * t238;
	t233 = t279 * t256;
	t226 = t255 * t236 + t258 * t241;
	t223 = -t258 * t236 + t255 * t241;
	t1 = [0, -t266 * qJD(4) - t278 * t267 - (t255 * t239 - t258 * t265) * t263 + (t255 * t274 - t258 * t262) * t277 + t264 * t226, 0, -t226, -t269 * t267 + ((-t251 * t276 + t252 * t266) * r_i_i_C(1) + (-t251 * t266 - t252 * t276) * r_i_i_C(2)) * qJD(5), 0; 0, t230 * qJD(4) + t278 * t268 - (-t258 * t239 - t255 * t265) * t263 + (-t255 * t262 - t258 * t274) * t277 + t264 * t223, 0, -t223, t269 * t268 + ((-t230 * t252 + t251 * t275) * r_i_i_C(1) + (t230 * t251 + t252 * t275) * r_i_i_C(2)) * qJD(5), 0; 0, t238 * qJD(4) + t278 * t233 - t264 * t234 + (-pkin(2) * t273 + t243 * t263) * t256, 0, t234, t269 * t233 + ((-t238 * t252 - t251 * t259) * r_i_i_C(1) + (t238 * t251 - t252 * t259) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:30
	% EndTime: 2019-10-09 21:24:31
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (575->89), mult. (1408->161), div. (0->0), fcn. (1519->13), ass. (0->58)
	t387 = pkin(12) + qJ(5);
	t385 = sin(t387);
	t386 = cos(t387);
	t394 = sin(qJ(6));
	t396 = cos(qJ(6));
	t403 = (t394 * r_i_i_C(1) + t396 * r_i_i_C(2)) * qJD(6);
	t408 = t396 * r_i_i_C(1) - t394 * r_i_i_C(2) + pkin(5);
	t425 = pkin(9) + r_i_i_C(3);
	t429 = (t408 * t385 - t425 * t386) * qJD(5) + t386 * t403;
	t392 = cos(pkin(6));
	t388 = sin(pkin(11));
	t395 = sin(qJ(2));
	t422 = cos(pkin(11));
	t424 = cos(qJ(2));
	t402 = t424 * t388 + t395 * t422;
	t373 = t402 * t392;
	t411 = t424 * t422;
	t417 = qJD(2) * t395;
	t427 = -qJD(2) * t411 + t388 * t417;
	t401 = -t395 * t388 + t411;
	t398 = t425 * t385 + t408 * t386 + cos(pkin(12)) * pkin(4) + pkin(3);
	t423 = pkin(2) * qJD(2);
	t389 = sin(pkin(10));
	t390 = sin(pkin(6));
	t421 = t389 * t390;
	t391 = cos(pkin(10));
	t420 = t390 * t391;
	t419 = t392 * t395;
	t416 = qJD(6) * t394;
	t415 = qJD(6) * t396;
	t370 = t427 * t392;
	t375 = t402 * qJD(2);
	t352 = t391 * t370 + t389 * t375;
	t354 = t389 * t370 - t391 * t375;
	t372 = t402 * t390;
	t363 = t372 * t386 + t392 * t385;
	t409 = -t372 * t385 + t392 * t386;
	t358 = t391 * t373 + t389 * t401;
	t359 = t389 * t373 - t391 * t401;
	t406 = -t358 * t385 - t386 * t420;
	t405 = -t358 * t386 + t385 * t420;
	t404 = t359 * t385 + t386 * t421;
	t349 = -t359 * t386 + t385 * t421;
	t400 = t401 * t392;
	t399 = qJD(2) * t373;
	t393 = -pkin(8) - qJ(4);
	t374 = t401 * qJD(2);
	t371 = t401 * t390;
	t369 = qJD(2) * t372;
	t368 = t427 * t390;
	t360 = -t389 * t400 - t391 * t402;
	t357 = -t389 * t402 + t391 * t400;
	t353 = -t391 * t374 + t389 * t399;
	t350 = -t389 * t374 - t391 * t399;
	t345 = t409 * qJD(5) - t368 * t386;
	t343 = t404 * qJD(5) + t354 * t386;
	t341 = t406 * qJD(5) - t352 * t386;
	t1 = [0, (t354 * t394 - t359 * t415) * r_i_i_C(1) + (t354 * t396 + t359 * t416) * r_i_i_C(2) - t354 * t393 - t359 * qJD(4) + (t389 * t419 - t424 * t391) * t423 + t398 * t353 - t429 * t360, 0, -t353, t425 * t343 - t404 * t403 + t408 * (-t349 * qJD(5) - t354 * t385), (-t343 * t394 - t353 * t396) * r_i_i_C(1) + (-t343 * t396 + t353 * t394) * r_i_i_C(2) + ((-t349 * t396 + t360 * t394) * r_i_i_C(1) + (t349 * t394 + t360 * t396) * r_i_i_C(2)) * qJD(6); 0, (-t352 * t394 + t358 * t415) * r_i_i_C(1) + (-t352 * t396 - t358 * t416) * r_i_i_C(2) + t352 * t393 + t358 * qJD(4) + (-t424 * t389 - t391 * t419) * t423 + t398 * t350 - t429 * t357, 0, -t350, t425 * t341 - t406 * t403 + t408 * (t405 * qJD(5) + t352 * t385), (-t341 * t394 - t350 * t396) * r_i_i_C(1) + (-t341 * t396 + t350 * t394) * r_i_i_C(2) + ((t357 * t394 + t396 * t405) * r_i_i_C(1) + (t357 * t396 - t394 * t405) * r_i_i_C(2)) * qJD(6); 0, (-t368 * t394 + t372 * t415) * r_i_i_C(1) + (-t368 * t396 - t372 * t416) * r_i_i_C(2) + t368 * t393 + t372 * qJD(4) - t390 * pkin(2) * t417 - t398 * t369 - t429 * t371, 0, t369, t425 * t345 - t409 * t403 + t408 * (-t363 * qJD(5) + t368 * t385), (-t345 * t394 + t369 * t396) * r_i_i_C(1) + (-t345 * t396 - t369 * t394) * r_i_i_C(2) + ((-t363 * t396 + t371 * t394) * r_i_i_C(1) + (t363 * t394 + t371 * t396) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end