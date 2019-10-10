% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:46
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR7_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR7_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR7_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:40
	% EndTime: 2019-10-10 09:46:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:41
	% EndTime: 2019-10-10 09:46:41
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (35->18), mult. (110->35), div. (0->0), fcn. (94->6), ass. (0->20)
	t136 = sin(pkin(6));
	t151 = t136 * (pkin(8) + r_i_i_C(3));
	t138 = sin(qJ(2));
	t139 = sin(qJ(1));
	t149 = t138 * t139;
	t141 = cos(qJ(1));
	t148 = t138 * t141;
	t140 = cos(qJ(2));
	t147 = t139 * t140;
	t146 = t140 * t141;
	t137 = cos(pkin(6));
	t145 = -t137 * t146 + t149;
	t144 = t137 * t147 + t148;
	t143 = t137 * t148 + t147;
	t142 = t137 * t149 - t146;
	t135 = t142 * qJD(1) + t145 * qJD(2);
	t134 = t144 * qJD(1) + t143 * qJD(2);
	t133 = t143 * qJD(1) + t144 * qJD(2);
	t132 = t145 * qJD(1) + t142 * qJD(2);
	t1 = [t135 * r_i_i_C(1) + t134 * r_i_i_C(2) + (-pkin(1) * t141 - t139 * t151) * qJD(1), t132 * r_i_i_C(1) + t133 * r_i_i_C(2), 0, 0, 0, 0; -t133 * r_i_i_C(1) + t132 * r_i_i_C(2) + (-pkin(1) * t139 + t141 * t151) * qJD(1), -t134 * r_i_i_C(1) + t135 * r_i_i_C(2), 0, 0, 0, 0; 0, (-r_i_i_C(1) * t138 - r_i_i_C(2) * t140) * t136 * qJD(2), 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:41
	% EndTime: 2019-10-10 09:46:41
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (84->29), mult. (244->45), div. (0->0), fcn. (222->6), ass. (0->24)
	t179 = sin(pkin(6));
	t198 = t179 * (pkin(8) + r_i_i_C(2));
	t197 = pkin(2) + r_i_i_C(1);
	t195 = r_i_i_C(3) + qJ(3);
	t181 = sin(qJ(2));
	t182 = sin(qJ(1));
	t194 = t182 * t181;
	t183 = cos(qJ(2));
	t193 = t182 * t183;
	t184 = cos(qJ(1));
	t192 = t184 * t181;
	t191 = t184 * t183;
	t190 = qJD(2) * t181;
	t180 = cos(pkin(6));
	t189 = t180 * t194;
	t188 = t180 * t191;
	t187 = qJD(2) * t180 + qJD(1);
	t186 = t180 * t193 + t192;
	t185 = t180 * t192 + t193;
	t174 = -qJD(1) * t189 - t182 * t190 + t187 * t191;
	t173 = t186 * qJD(1) + t185 * qJD(2);
	t172 = t185 * qJD(1) + t186 * qJD(2);
	t171 = -qJD(1) * t188 - qJD(2) * t191 + t187 * t194;
	t1 = [-(-t188 + t194) * qJD(3) - t197 * t174 - t195 * t173 + (-t184 * pkin(1) - t182 * t198) * qJD(1), -(t189 - t191) * qJD(3) - t195 * t172 + t197 * t171, -t171, 0, 0, 0; t186 * qJD(3) - t197 * t172 - t195 * t171 + (-t182 * pkin(1) + t184 * t198) * qJD(1), t185 * qJD(3) - t197 * t173 + t195 * t174, t173, 0, 0, 0; 0, (t181 * qJD(3) + (-t197 * t181 + t195 * t183) * qJD(2)) * t179, t179 * t190, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:41
	% EndTime: 2019-10-10 09:46:41
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (107->34), mult. (305->51), div. (0->0), fcn. (276->6), ass. (0->26)
	t148 = sin(pkin(6));
	t169 = t148 * (pkin(8) - r_i_i_C(3) - qJ(4));
	t168 = r_i_i_C(1) + qJ(3);
	t150 = sin(qJ(2));
	t151 = sin(qJ(1));
	t167 = t151 * t150;
	t152 = cos(qJ(2));
	t166 = t151 * t152;
	t153 = cos(qJ(1));
	t165 = t153 * t150;
	t164 = t153 * t152;
	t163 = qJD(1) * t148;
	t162 = qJD(2) * t150;
	t161 = t148 * qJD(4);
	t160 = pkin(2) + pkin(3) - r_i_i_C(2);
	t149 = cos(pkin(6));
	t158 = t149 * t167;
	t157 = t149 * t164;
	t156 = qJD(2) * t149 + qJD(1);
	t155 = t149 * t166 + t165;
	t154 = t149 * t165 + t166;
	t143 = -qJD(1) * t158 - t151 * t162 + t156 * t164;
	t142 = t155 * qJD(1) + t154 * qJD(2);
	t141 = t154 * qJD(1) + t155 * qJD(2);
	t140 = -qJD(1) * t157 - qJD(2) * t164 + t156 * t167;
	t1 = [-t153 * t161 - (-t157 + t167) * qJD(3) - t168 * t142 - t160 * t143 + (-t153 * pkin(1) - t151 * t169) * qJD(1), -(t158 - t164) * qJD(3) - t168 * t141 + t160 * t140, -t140, -t153 * t163, 0, 0; -t151 * t161 + t155 * qJD(3) - t168 * t140 - t160 * t141 + (-t151 * pkin(1) + t153 * t169) * qJD(1), t154 * qJD(3) - t160 * t142 + t168 * t143, t142, -t151 * t163, 0, 0; 0, (t150 * qJD(3) + (-t160 * t150 + t168 * t152) * qJD(2)) * t148, t148 * t162, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:41
	% EndTime: 2019-10-10 09:46:42
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (217->60), mult. (636->95), div. (0->0), fcn. (594->8), ass. (0->43)
	t241 = sin(qJ(5));
	t244 = cos(qJ(5));
	t271 = r_i_i_C(2) * t244;
	t250 = r_i_i_C(1) * t241 + t271;
	t272 = t250 * qJD(5) - qJD(3);
	t270 = -pkin(4) - qJ(3);
	t269 = pkin(8) - qJ(4);
	t239 = sin(pkin(6));
	t243 = sin(qJ(1));
	t268 = t239 * t243;
	t245 = cos(qJ(2));
	t267 = t239 * t245;
	t246 = cos(qJ(1));
	t266 = t239 * t246;
	t242 = sin(qJ(2));
	t265 = t243 * t242;
	t264 = t243 * t245;
	t263 = t246 * t242;
	t262 = t246 * t245;
	t261 = qJD(1) * t243;
	t260 = qJD(1) * t246;
	t259 = qJD(2) * t242;
	t258 = -r_i_i_C(3) - pkin(9) - pkin(3) - pkin(2);
	t240 = cos(pkin(6));
	t257 = t240 * t265;
	t256 = t240 * t262;
	t255 = t239 * t261;
	t254 = t239 * t260;
	t253 = t239 * t259;
	t252 = qJD(2) * t240 + qJD(1);
	t251 = -t244 * r_i_i_C(1) + t241 * r_i_i_C(2);
	t231 = t240 * t264 + t263;
	t249 = t240 * t263 + t264;
	t248 = -t251 - t270;
	t233 = t241 * t255;
	t229 = -t256 + t265;
	t228 = -qJD(1) * t257 - t243 * t259 + t252 * t262;
	t227 = t231 * qJD(1) + t249 * qJD(2);
	t226 = t249 * qJD(1) + t231 * qJD(2);
	t225 = -qJD(1) * t256 - qJD(2) * t262 + t252 * t265;
	t224 = -t241 * t254 - t225 * t244 + (-t231 * t241 - t244 * t268) * qJD(5);
	t223 = -t244 * t254 + t225 * t241 + (-t231 * t244 + t241 * t268) * qJD(5);
	t1 = [-pkin(1) * t260 + t233 * r_i_i_C(1) + t272 * t229 - t248 * t227 + ((t251 * qJD(5) - qJD(4)) * t246 + (-t269 + t271) * t261) * t239 + t258 * t228, -t272 * (-t257 + t262) - t248 * t226 - t258 * t225, -t225, -t254, t223 * r_i_i_C(1) - t224 * r_i_i_C(2), 0; -qJD(4) * t268 + t224 * r_i_i_C(1) + t223 * r_i_i_C(2) + t231 * qJD(3) + t270 * t225 + (-pkin(1) * t243 + t269 * t266) * qJD(1) + t258 * t226, t258 * t227 + t248 * t228 - t249 * t272, t227, -t255, (-t227 * t241 - t244 * t255) * r_i_i_C(1) + (-t227 * t244 + t233) * r_i_i_C(2) + ((-t229 * t244 - t241 * t266) * r_i_i_C(1) + (t229 * t241 - t244 * t266) * r_i_i_C(2)) * qJD(5), 0; 0, (-t272 * t242 + (t258 * t242 + t248 * t245) * qJD(2)) * t239, t253, 0, -t250 * t253 + ((t240 * t241 + t244 * t267) * r_i_i_C(1) + (t240 * t244 - t241 * t267) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:46:43
	% EndTime: 2019-10-10 09:46:43
	% DurationCPUTime: 0.66s
	% Computational Cost: add. (516->98), mult. (1526->155), div. (0->0), fcn. (1516->10), ass. (0->66)
	t363 = cos(pkin(6));
	t370 = cos(qJ(2));
	t371 = cos(qJ(1));
	t399 = t370 * t371;
	t391 = t363 * t399;
	t366 = sin(qJ(2));
	t367 = sin(qJ(1));
	t402 = t366 * t367;
	t351 = -t391 + t402;
	t365 = sin(qJ(5));
	t369 = cos(qJ(5));
	t362 = sin(pkin(6));
	t403 = t362 * t371;
	t345 = t351 * t369 + t365 * t403;
	t400 = t367 * t370;
	t401 = t366 * t371;
	t352 = t363 * t401 + t400;
	t364 = sin(qJ(6));
	t368 = cos(qJ(6));
	t413 = t345 * t364 - t352 * t368;
	t412 = t345 * t368 + t352 * t364;
	t384 = r_i_i_C(1) * t368 - r_i_i_C(2) * t364;
	t377 = t384 * qJD(6);
	t353 = t363 * t400 + t401;
	t340 = t353 * qJD(1) + t352 * qJD(2);
	t381 = -t351 * t365 + t369 * t403;
	t398 = qJD(1) * t362;
	t390 = t367 * t398;
	t335 = t381 * qJD(5) + t340 * t369 - t365 * t390;
	t382 = pkin(5) + t384;
	t410 = -pkin(4) - qJ(3);
	t411 = pkin(10) + r_i_i_C(3);
	t374 = t411 * t365 + t382 * t369 - t410;
	t409 = pkin(8) - qJ(4);
	t406 = t362 * t366;
	t405 = t362 * t367;
	t404 = t362 * t370;
	t397 = qJD(2) * t366;
	t396 = qJD(2) * t370;
	t395 = t362 * qJD(4);
	t394 = -pkin(2) - pkin(3) - pkin(9);
	t393 = t365 * t405;
	t392 = t363 * t402;
	t389 = t371 * t398;
	t388 = t362 * t396;
	t387 = t362 * t397;
	t385 = qJD(2) * t363 + qJD(1);
	t383 = -t364 * r_i_i_C(1) - t368 * r_i_i_C(2);
	t380 = -t353 * t365 - t369 * t405;
	t379 = t363 * t365 + t369 * t404;
	t378 = -t363 * t369 + t365 * t404;
	t376 = qJD(6) * t383;
	t375 = -t383 - t394;
	t373 = -t345 * qJD(5) - t340 * t365 - t369 * t390;
	t372 = qJD(3) + t369 * t376 + (-t382 * t365 + t411 * t369) * qJD(5);
	t354 = -t392 + t399;
	t348 = t353 * t369 - t393;
	t343 = t378 * qJD(5) + t369 * t387;
	t341 = -qJD(1) * t392 - t367 * t397 + t385 * t399;
	t339 = t352 * qJD(1) + t353 * qJD(2);
	t338 = -qJD(1) * t391 - t371 * t396 + t385 * t402;
	t333 = t380 * qJD(5) - t338 * t369 - t365 * t389;
	t332 = -t338 * t365 - qJD(5) * t393 + (qJD(5) * t353 + t389) * t369;
	t331 = t333 * t368 - t339 * t364 + (-t348 * t364 + t354 * t368) * qJD(6);
	t330 = -t333 * t364 - t339 * t368 + (-t348 * t368 - t354 * t364) * qJD(6);
	t1 = [-t371 * t395 - t351 * qJD(3) + t410 * t340 + t411 * t373 - t382 * t335 + (t413 * r_i_i_C(1) + t412 * r_i_i_C(2)) * qJD(6) - t375 * t341 + (-t371 * pkin(1) - t409 * t405) * qJD(1), t375 * t338 - t374 * t339 - t353 * t377 + t372 * t354, -t338, -t389, -t382 * t332 + t411 * t333 + t380 * t376, r_i_i_C(1) * t330 - r_i_i_C(2) * t331; -t367 * t395 + t333 * pkin(5) + t331 * r_i_i_C(1) + t330 * r_i_i_C(2) + t353 * qJD(3) + t410 * t338 + t411 * t332 + t394 * t339 + (-pkin(1) * t367 + t409 * t403) * qJD(1), -t375 * t340 + t374 * t341 - t351 * t377 + t372 * t352, t340, -t390, t411 * t335 + t382 * t373 + t381 * t376, (-t335 * t364 + t341 * t368) * r_i_i_C(1) + (-t335 * t368 - t341 * t364) * r_i_i_C(2) + (-t412 * r_i_i_C(1) + t413 * r_i_i_C(2)) * qJD(6); 0, ((t374 * qJD(2) + t377) * t370 + (-t375 * qJD(2) + t372) * t366) * t362, t387, 0, t411 * t343 + t378 * t376 + t382 * (t379 * qJD(5) - t365 * t387), (-t343 * t364 + t368 * t388) * r_i_i_C(1) + (-t343 * t368 - t364 * t388) * r_i_i_C(2) + ((-t364 * t406 + t368 * t379) * r_i_i_C(1) + (-t364 * t379 - t368 * t406) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end