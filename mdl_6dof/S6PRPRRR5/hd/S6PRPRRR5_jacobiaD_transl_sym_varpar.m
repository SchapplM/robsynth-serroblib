% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S6PRPRRR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% JaD_transl [3x6]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:01
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR5_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_transl_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_transl_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR5_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR5_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobiaD_transl_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:16
	% EndTime: 2019-10-09 22:01:16
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
	% StartTime: 2019-10-09 22:01:17
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (24->13), mult. (82->26), div. (0->0), fcn. (72->6), ass. (0->16)
	t135 = -pkin(2) + r_i_i_C(2);
	t134 = r_i_i_C(3) + qJ(3);
	t126 = cos(pkin(6));
	t127 = sin(qJ(2));
	t133 = t126 * t127;
	t128 = cos(qJ(2));
	t132 = t126 * t128;
	t131 = qJD(2) * t134;
	t123 = sin(pkin(11));
	t125 = cos(pkin(11));
	t130 = t123 * t128 + t125 * t133;
	t129 = t123 * t133 - t125 * t128;
	t124 = sin(pkin(6));
	t122 = t129 * qJD(2);
	t120 = t130 * qJD(2);
	t1 = [0, -t129 * qJD(3) - t135 * t122 - (t123 * t132 + t125 * t127) * t131, -t122, 0, 0, 0; 0, t130 * qJD(3) + t135 * t120 - (t123 * t127 - t125 * t132) * t131, t120, 0, 0, 0; 0, (qJD(3) * t127 + (t135 * t127 + t134 * t128) * qJD(2)) * t124, t124 * qJD(2) * t127, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:17
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (67->29), mult. (228->61), div. (0->0), fcn. (208->8), ass. (0->28)
	t178 = sin(pkin(6));
	t181 = sin(qJ(4));
	t199 = t178 * t181;
	t183 = cos(qJ(4));
	t198 = t178 * t183;
	t184 = cos(qJ(2));
	t197 = t178 * t184;
	t179 = cos(pkin(11));
	t196 = t179 * t184;
	t180 = cos(pkin(6));
	t182 = sin(qJ(2));
	t195 = t180 * t182;
	t194 = t180 * t184;
	t193 = qJD(2) * t182;
	t192 = -pkin(2) - pkin(8) - r_i_i_C(3);
	t177 = sin(pkin(11));
	t191 = t177 * t193;
	t190 = t178 * t193;
	t189 = qJD(2) * t196;
	t188 = t183 * r_i_i_C(1) - t181 * r_i_i_C(2);
	t187 = t181 * r_i_i_C(1) + t183 * r_i_i_C(2) + qJ(3);
	t186 = t177 * t184 + t179 * t195;
	t173 = t177 * t194 + t179 * t182;
	t185 = t188 * qJD(4) + qJD(3);
	t171 = t177 * t182 - t179 * t194;
	t170 = -t180 * t191 + t189;
	t168 = t186 * qJD(2);
	t1 = [0, t185 * (-t177 * t195 + t196) + t192 * t170 - t187 * t173 * qJD(2), t170, t188 * t170 + ((-t173 * t181 - t177 * t198) * r_i_i_C(1) + (-t173 * t183 + t177 * t199) * r_i_i_C(2)) * qJD(4), 0, 0; 0, t185 * t186 + t192 * t168 - t187 * (-t180 * t189 + t191), t168, t188 * t168 + ((-t171 * t181 + t179 * t198) * r_i_i_C(1) + (-t171 * t183 - t179 * t199) * r_i_i_C(2)) * qJD(4), 0, 0; 0, (t185 * t182 + (t192 * t182 + t187 * t184) * qJD(2)) * t178, t190, t188 * t190 + ((-t180 * t183 + t181 * t197) * r_i_i_C(1) + (t180 * t181 + t183 * t197) * r_i_i_C(2)) * qJD(4), 0, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:17
	% EndTime: 2019-10-09 22:01:17
	% DurationCPUTime: 0.17s
	% Computational Cost: add. (194->42), mult. (385->79), div. (0->0), fcn. (352->10), ass. (0->40)
	t221 = qJ(4) + qJ(5);
	t218 = sin(t221);
	t220 = qJD(4) + qJD(5);
	t253 = t218 * t220;
	t219 = cos(t221);
	t252 = t219 * t220;
	t223 = sin(pkin(6));
	t251 = t220 * t223;
	t228 = cos(qJ(4));
	t250 = t223 * t228;
	t229 = cos(qJ(2));
	t249 = t223 * t229;
	t224 = cos(pkin(11));
	t248 = t224 * t229;
	t225 = cos(pkin(6));
	t227 = sin(qJ(2));
	t247 = t225 * t227;
	t246 = t225 * t229;
	t222 = sin(pkin(11));
	t214 = t222 * t246 + t224 * t227;
	t237 = qJD(2) * t248;
	t242 = qJD(2) * t227;
	t239 = t222 * t242;
	t211 = -t225 * t239 + t237;
	t235 = t222 * t251 - t211;
	t245 = (-t214 * t253 - t235 * t219) * r_i_i_C(1) + (-t214 * t252 + t235 * t218) * r_i_i_C(2);
	t212 = t222 * t227 - t224 * t246;
	t234 = t222 * t229 + t224 * t247;
	t209 = t234 * qJD(2);
	t236 = t224 * t251 + t209;
	t244 = (-t212 * t253 + t236 * t219) * r_i_i_C(1) + (-t212 * t252 - t236 * t218) * r_i_i_C(2);
	t238 = t223 * t242;
	t233 = -t220 * t225 + t238;
	t240 = t220 * t249;
	t243 = (t218 * t240 + t233 * t219) * r_i_i_C(1) + (-t233 * t218 + t219 * t240) * r_i_i_C(2);
	t241 = -pkin(2) - r_i_i_C(3) - pkin(9) - pkin(8);
	t226 = sin(qJ(4));
	t232 = pkin(4) * t226 + r_i_i_C(1) * t218 + r_i_i_C(2) * t219 + qJ(3);
	t231 = pkin(4) * qJD(4) * t228 + qJD(3) + (r_i_i_C(1) * t219 - r_i_i_C(2) * t218) * t220;
	t1 = [0, t241 * t211 + t231 * (-t222 * t247 + t248) - t232 * t214 * qJD(2), t211, (t211 * t228 + (-t214 * t226 - t222 * t250) * qJD(4)) * pkin(4) + t245, t245, 0; 0, t241 * t209 + t231 * t234 - t232 * (-t225 * t237 + t239), t209, (t209 * t228 + (-t212 * t226 + t224 * t250) * qJD(4)) * pkin(4) + t244, t244, 0; 0, (t231 * t227 + (t241 * t227 + t232 * t229) * qJD(2)) * t223, t238, (t228 * t238 + (-t225 * t228 + t226 * t249) * qJD(4)) * pkin(4) + t243, t243, 0;];
	JaD_transl = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_transl_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:01:18
	% EndTime: 2019-10-09 22:01:19
	% DurationCPUTime: 0.79s
	% Computational Cost: add. (607->81), mult. (1090->146), div. (0->0), fcn. (1076->12), ass. (0->60)
	t358 = sin(qJ(6));
	t361 = cos(qJ(6));
	t405 = t358 * r_i_i_C(1) + t361 * r_i_i_C(2);
	t402 = t405 * qJD(6);
	t399 = pkin(10) + r_i_i_C(3);
	t376 = t361 * r_i_i_C(1) - t358 * r_i_i_C(2);
	t403 = pkin(5) + t376;
	t371 = t376 * qJD(6);
	t353 = qJ(4) + qJ(5);
	t350 = sin(t353);
	t351 = cos(t353);
	t359 = sin(qJ(4));
	t400 = t359 * pkin(4) + t350 * t403 - t399 * t351 + qJ(3);
	t352 = qJD(4) + qJD(5);
	t396 = t350 * t352;
	t395 = t351 * t352;
	t354 = sin(pkin(11));
	t355 = sin(pkin(6));
	t394 = t354 * t355;
	t356 = cos(pkin(11));
	t393 = t355 * t356;
	t360 = sin(qJ(2));
	t392 = t355 * t360;
	t362 = cos(qJ(4));
	t391 = t355 * t362;
	t363 = cos(qJ(2));
	t390 = t355 * t363;
	t357 = cos(pkin(6));
	t389 = t357 * t360;
	t388 = t357 * t363;
	t387 = qJD(2) * t360;
	t386 = qJD(2) * t363;
	t383 = t350 * t394;
	t382 = t350 * t390;
	t381 = t355 * t386;
	t380 = t354 * t387;
	t379 = t355 * t387;
	t378 = t356 * t386;
	t343 = t354 * t363 + t356 * t389;
	t339 = t343 * qJD(2);
	t374 = t352 * t393 + t339;
	t372 = t357 * t350 + t351 * t390;
	t344 = t354 * t388 + t356 * t360;
	t370 = -pkin(2) - pkin(9) - pkin(8) - t405;
	t341 = -t357 * t380 + t378;
	t321 = -t341 * t350 - t344 * t395 + t352 * t383;
	t368 = -t402 * (t344 * t351 - t383) + t403 * (-t344 * t396 + (-t352 * t394 + t341) * t351) - t399 * t321;
	t342 = t354 * t360 - t356 * t388;
	t323 = t342 * t395 + t374 * t350;
	t367 = -t402 * (t342 * t351 + t350 * t393) + t403 * (-t342 * t396 + t374 * t351) + t399 * t323;
	t328 = -t350 * t379 + t372 * t352;
	t366 = t402 * t372 + t403 * (t352 * t382 + (-t352 * t357 + t379) * t351) - t399 * t328;
	t365 = qJD(4) * t362 * pkin(4) + qJD(3) - t350 * t402 + (t399 * t350 + t351 * t403) * t352;
	t345 = -t354 * t389 + t356 * t363;
	t340 = t344 * qJD(2);
	t338 = -t357 * t378 + t380;
	t337 = t357 * t351 - t382;
	t333 = t342 * t350 - t351 * t393;
	t331 = t344 * t350 + t351 * t394;
	t1 = [0, -t340 * t400 + t370 * t341 - t344 * t371 + t365 * t345, t341, (t341 * t362 + (-t344 * t359 - t354 * t391) * qJD(4)) * pkin(4) + t368, t368, (t321 * t358 - t340 * t361) * r_i_i_C(1) + (t321 * t361 + t340 * t358) * r_i_i_C(2) + ((-t331 * t361 - t345 * t358) * r_i_i_C(1) + (t331 * t358 - t345 * t361) * r_i_i_C(2)) * qJD(6); 0, -t338 * t400 + t370 * t339 - t342 * t371 + t365 * t343, t339, (t339 * t362 + (-t342 * t359 + t356 * t391) * qJD(4)) * pkin(4) + t367, t367, (-t323 * t358 - t338 * t361) * r_i_i_C(1) + (-t323 * t361 + t338 * t358) * r_i_i_C(2) + ((-t333 * t361 - t343 * t358) * r_i_i_C(1) + (t333 * t358 - t343 * t361) * r_i_i_C(2)) * qJD(6); 0, ((t400 * qJD(2) + t371) * t363 + (t370 * qJD(2) + t365) * t360) * t355, t379, (t362 * t379 + (-t357 * t362 + t359 * t390) * qJD(4)) * pkin(4) + t366, t366, (t328 * t358 + t361 * t381) * r_i_i_C(1) + (t328 * t361 - t358 * t381) * r_i_i_C(2) + ((-t337 * t361 - t358 * t392) * r_i_i_C(1) + (t337 * t358 - t361 * t392) * r_i_i_C(2)) * qJD(6);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,6);
end