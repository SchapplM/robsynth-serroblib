% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6RRPPRR5
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
% Datum: 2019-10-10 09:43
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiaD_transl_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:06
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
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:07
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
	% StartTime: 2019-10-10 09:43:07
	% EndTime: 2019-10-10 09:43:07
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-10 09:43:06
	% EndTime: 2019-10-10 09:43:07
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (107->34), mult. (305->51), div. (0->0), fcn. (276->6), ass. (0->26)
	t147 = sin(pkin(6));
	t168 = t147 * (pkin(8) - r_i_i_C(3) - qJ(4));
	t167 = r_i_i_C(2) + qJ(3);
	t149 = sin(qJ(2));
	t150 = sin(qJ(1));
	t166 = t150 * t149;
	t151 = cos(qJ(2));
	t165 = t150 * t151;
	t152 = cos(qJ(1));
	t164 = t152 * t149;
	t163 = t152 * t151;
	t162 = qJD(1) * t147;
	t161 = qJD(2) * t149;
	t160 = t147 * qJD(4);
	t159 = pkin(2) + pkin(3) + r_i_i_C(1);
	t148 = cos(pkin(6));
	t157 = t148 * t166;
	t156 = t148 * t163;
	t155 = qJD(2) * t148 + qJD(1);
	t154 = t148 * t165 + t164;
	t153 = t148 * t164 + t165;
	t142 = -qJD(1) * t157 - t150 * t161 + t155 * t163;
	t141 = t154 * qJD(1) + t153 * qJD(2);
	t140 = t153 * qJD(1) + t154 * qJD(2);
	t139 = -qJD(1) * t156 - qJD(2) * t163 + t155 * t166;
	t1 = [-t152 * t160 - (-t156 + t166) * qJD(3) - t167 * t141 - t159 * t142 + (-t152 * pkin(1) - t150 * t168) * qJD(1), -(t157 - t163) * qJD(3) - t167 * t140 + t159 * t139, -t139, -t152 * t162, 0, 0; -t150 * t160 + t154 * qJD(3) - t167 * t139 - t159 * t140 + (-t150 * pkin(1) + t152 * t168) * qJD(1), t153 * qJD(3) - t159 * t141 + t167 * t142, t141, -t150 * t162, 0, 0; 0, (t149 * qJD(3) + (-t159 * t149 + t167 * t151) * qJD(2)) * t147, t147 * t161, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:07
	% EndTime: 2019-10-10 09:43:07
	% DurationCPUTime: 0.31s
	% Computational Cost: add. (217->64), mult. (636->104), div. (0->0), fcn. (594->8), ass. (0->45)
	t244 = sin(qJ(5));
	t247 = cos(qJ(5));
	t275 = t247 * r_i_i_C(2);
	t253 = t244 * r_i_i_C(1) + t275;
	t251 = qJD(5) * t253;
	t274 = pkin(8) - qJ(4);
	t242 = sin(pkin(6));
	t246 = sin(qJ(1));
	t273 = t242 * t246;
	t272 = t242 * t247;
	t249 = cos(qJ(1));
	t271 = t242 * t249;
	t245 = sin(qJ(2));
	t270 = t245 * t246;
	t269 = t245 * t249;
	t248 = cos(qJ(2));
	t268 = t246 * t248;
	t267 = t248 * t249;
	t266 = qJD(1) * t246;
	t265 = qJD(1) * t249;
	t264 = qJD(2) * t245;
	t263 = qJD(2) * t248;
	t243 = cos(pkin(6));
	t233 = t243 * t269 + t268;
	t262 = qJD(5) * t233;
	t261 = -pkin(2) - pkin(3) - pkin(4);
	t260 = r_i_i_C(3) + pkin(9) - qJ(3);
	t259 = t243 * t270;
	t258 = t243 * t267;
	t257 = t242 * t266;
	t256 = t242 * t265;
	t255 = qJD(2) * t243 + qJD(1);
	t254 = -t247 * r_i_i_C(1) + t244 * r_i_i_C(2);
	t252 = t243 * t268 + t269;
	t250 = -t254 - t261;
	t236 = t244 * t257;
	t235 = -t259 + t267;
	t232 = -t258 + t270;
	t231 = -qJD(1) * t259 - t246 * t264 + t255 * t267;
	t230 = t252 * qJD(1) + t233 * qJD(2);
	t229 = t233 * qJD(1) + t252 * qJD(2);
	t228 = -qJD(1) * t258 - t249 * t263 + t255 * t270;
	t227 = -t244 * t256 - t229 * t247 + (-t235 * t244 - t246 * t272) * qJD(5);
	t226 = -t247 * t256 + t229 * t244 + (-t235 * t247 + t244 * t273) * qJD(5);
	t1 = [(t244 * t262 + t236) * r_i_i_C(1) + t262 * t275 - t232 * qJD(3) - pkin(1) * t265 - t250 * t231 + t260 * t230 + ((t254 * qJD(5) - qJD(4)) * t249 + (-t274 + t275) * t266) * t242, qJD(3) * t235 + t250 * t228 + t260 * t229 + t251 * t252, -t228, -t256, r_i_i_C(1) * t226 - r_i_i_C(2) * t227, 0; -qJD(4) * t273 + t227 * r_i_i_C(1) + t226 * r_i_i_C(2) + t252 * qJD(3) + t261 * t229 + t260 * t228 + (-pkin(1) * t246 + t274 * t271) * qJD(1), qJD(3) * t233 - t250 * t230 - t260 * t231 + t232 * t251, t230, -t257, (-t231 * t244 - t247 * t257) * r_i_i_C(1) + (-t231 * t247 + t236) * r_i_i_C(2) + ((-t233 * t247 - t244 * t271) * r_i_i_C(1) + (t233 * t244 - t247 * t271) * r_i_i_C(2)) * qJD(5), 0; 0, (qJD(3) * t245 - t248 * t251 + (-t250 * t245 - t260 * t248) * qJD(2)) * t242, t242 * t264, 0, -t253 * t242 * t263 + ((t243 * t244 - t245 * t272) * r_i_i_C(1) + (t242 * t244 * t245 + t243 * t247) * r_i_i_C(2)) * qJD(5), 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:43:09
	% EndTime: 2019-10-10 09:43:09
	% DurationCPUTime: 0.73s
	% Computational Cost: add. (516->98), mult. (1526->159), div. (0->0), fcn. (1516->10), ass. (0->66)
	t365 = sin(qJ(2));
	t366 = sin(qJ(1));
	t369 = cos(qJ(2));
	t370 = cos(qJ(1));
	t408 = cos(pkin(6));
	t388 = t370 * t408;
	t351 = t365 * t388 + t366 * t369;
	t364 = sin(qJ(5));
	t368 = cos(qJ(5));
	t362 = sin(pkin(6));
	t402 = t362 * t370;
	t344 = t351 * t368 + t364 * t402;
	t385 = t369 * t388;
	t401 = t366 * t365;
	t350 = -t385 + t401;
	t363 = sin(qJ(6));
	t367 = cos(qJ(6));
	t417 = t344 * t363 + t350 * t367;
	t416 = t344 * t367 - t350 * t363;
	t383 = t363 * r_i_i_C(1) + t367 * r_i_i_C(2);
	t378 = qJD(6) * t383;
	t384 = -t367 * r_i_i_C(1) + t363 * r_i_i_C(2);
	t381 = pkin(5) - t384;
	t411 = pkin(10) + r_i_i_C(3);
	t415 = (t381 * t364 - t411 * t368) * qJD(5) + t368 * t378;
	t382 = qJD(2) * t408 + qJD(1);
	t389 = t366 * t408;
	t386 = t365 * t389;
	t398 = qJD(2) * t365;
	t400 = t370 * t369;
	t340 = -qJD(1) * t386 - t366 * t398 + t382 * t400;
	t380 = -t351 * t364 + t368 * t402;
	t399 = qJD(1) * t362;
	t393 = t366 * t399;
	t334 = t380 * qJD(5) + t340 * t368 - t364 * t393;
	t395 = -pkin(2) - pkin(3) - pkin(4);
	t412 = t411 * t364 + t381 * t368 - t395;
	t410 = pkin(8) - qJ(4);
	t409 = pkin(9) - qJ(3);
	t405 = t362 * t366;
	t404 = t362 * t368;
	t403 = t362 * t369;
	t397 = qJD(2) * t369;
	t396 = t362 * qJD(4);
	t394 = t364 * t405;
	t392 = t370 * t399;
	t391 = t362 * t397;
	t390 = t362 * t398;
	t353 = -t386 + t400;
	t379 = -t353 * t364 - t366 * t404;
	t377 = t383 + t409;
	t376 = -t362 * t365 * t364 - t408 * t368;
	t349 = -t408 * t364 + t365 * t404;
	t375 = t370 * t365 + t369 * t389;
	t374 = t384 * qJD(6) + qJD(3);
	t372 = -qJD(5) * t344 - t340 * t364 - t368 * t393;
	t347 = t353 * t368 - t394;
	t342 = qJD(5) * t376 + t368 * t391;
	t339 = qJD(1) * t375 + qJD(2) * t351;
	t338 = qJD(1) * t351 + qJD(2) * t375;
	t337 = -qJD(1) * t385 - t370 * t397 + t382 * t401;
	t332 = qJD(5) * t379 - t338 * t368 - t364 * t392;
	t331 = -t338 * t364 - qJD(5) * t394 + (qJD(5) * t353 + t392) * t368;
	t330 = t332 * t367 + t337 * t363 + (-t347 * t363 - t367 * t375) * qJD(6);
	t329 = -t332 * t363 + t337 * t367 + (-t347 * t367 + t363 * t375) * qJD(6);
	t1 = [-t370 * t396 - t350 * qJD(3) - t381 * t334 + t377 * t339 + t411 * t372 + (t417 * r_i_i_C(1) + t416 * r_i_i_C(2)) * qJD(6) + t395 * t340 + (-t370 * pkin(1) - t410 * t405) * qJD(1), t412 * t337 + t377 * t338 + t374 * t353 + t375 * t415, -t337, -t392, -t381 * t331 + t411 * t332 - t379 * t378, t329 * r_i_i_C(1) - t330 * r_i_i_C(2); -t366 * t396 + t332 * pkin(5) + t330 * r_i_i_C(1) + t329 * r_i_i_C(2) + t375 * qJD(3) + t409 * t337 + t411 * t331 + t395 * t338 + (-pkin(1) * t366 + t410 * t402) * qJD(1), -t339 * t412 - t377 * t340 + t415 * t350 + t374 * t351, t339, -t393, t411 * t334 + t381 * t372 - t380 * t378, (-t334 * t363 - t339 * t367) * r_i_i_C(1) + (-t334 * t367 + t339 * t363) * r_i_i_C(2) + (-t416 * r_i_i_C(1) + t417 * r_i_i_C(2)) * qJD(6); 0, (t374 * t365 - t415 * t369 + (-t365 * t412 - t369 * t377) * qJD(2)) * t362, t390, 0, t411 * t342 - t376 * t378 + t381 * (-qJD(5) * t349 - t364 * t391), (-t342 * t363 - t367 * t390) * r_i_i_C(1) + (-t342 * t367 + t363 * t390) * r_i_i_C(2) + ((-t349 * t367 - t363 * t403) * r_i_i_C(1) + (t349 * t363 - t367 * t403) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end