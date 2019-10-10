% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PRRPRR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiaD_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(12));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(12)) * cos(pkin(6));
	t37 = t39 * t42 - t41 * t45;
	t34 = 0.1e1 / t37 ^ 2;
	t49 = qJD(2) * t34;
	t36 = t39 * t41 + t42 * t45;
	t33 = t36 ^ 2;
	t30 = t33 * t34 + 0.1e1;
	t46 = t37 * t49;
	t47 = t36 / t37 * t49;
	t48 = (t33 * t47 + t36 * t46) / t30 ^ 2;
	t28 = 0.1e1 / t30;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:27
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (1326->63), mult. (4214->161), div. (275->12), fcn. (5416->13), ass. (0->79)
	t164 = sin(pkin(7));
	t162 = t164 ^ 2;
	t206 = 0.2e1 * t162;
	t166 = cos(pkin(12));
	t163 = sin(pkin(12));
	t170 = sin(qJ(2));
	t168 = cos(pkin(6));
	t172 = cos(qJ(2));
	t190 = t168 * t172;
	t182 = -t163 * t170 + t166 * t190;
	t165 = sin(pkin(6));
	t167 = cos(pkin(7));
	t195 = t165 * t167;
	t146 = t182 * t164 + t166 * t195;
	t196 = t164 * t172;
	t156 = -t165 * t196 + t168 * t167;
	t141 = atan2(t146, t156);
	t136 = sin(t141);
	t137 = cos(t141);
	t123 = t136 * t146 + t137 * t156;
	t120 = 0.1e1 / t123;
	t169 = sin(qJ(3));
	t171 = cos(qJ(3));
	t191 = t168 * t170;
	t180 = t163 * t191 - t166 * t172;
	t181 = t163 * t190 + t166 * t170;
	t187 = t163 * t164 * t165;
	t183 = -t167 * t181 + t187;
	t135 = t183 * t169 - t171 * t180;
	t131 = 0.1e1 / t135;
	t153 = 0.1e1 / t156;
	t121 = 0.1e1 / t123 ^ 2;
	t132 = 0.1e1 / t135 ^ 2;
	t154 = 0.1e1 / t156 ^ 2;
	t147 = t163 * t195 + t164 * t181;
	t145 = t147 ^ 2;
	t119 = t121 * t145 + 0.1e1;
	t152 = t180 * qJD(2);
	t157 = -t163 * t172 - t166 * t191;
	t150 = t157 * qJD(2);
	t194 = t165 * t170;
	t198 = t146 * t154;
	t185 = t194 * t198;
	t144 = t146 ^ 2;
	t140 = t144 * t154 + 0.1e1;
	t138 = 0.1e1 / t140;
	t199 = t138 * t164;
	t115 = (-qJD(2) * t185 + t150 * t153) * t199;
	t184 = -t136 * t156 + t137 * t146;
	t189 = qJD(2) * t165;
	t186 = t170 * t189;
	t112 = (t136 * t150 + t137 * t186) * t164 + t184 * t115;
	t204 = t112 * t120 * t121;
	t205 = (-t121 * t147 * t152 * t164 - t145 * t204) / t119 ^ 2;
	t192 = t167 * t171;
	t197 = t180 * t169;
	t134 = -t171 * t187 + t181 * t192 - t197;
	t130 = t134 ^ 2;
	t127 = t130 * t132 + 0.1e1;
	t151 = t181 * qJD(2);
	t193 = t167 * t169;
	t129 = t152 * t193 - t151 * t171 + (t183 * t171 + t197) * qJD(3);
	t201 = t129 * t131 * t132;
	t128 = t135 * qJD(3) - t151 * t169 - t152 * t192;
	t202 = t128 * t132;
	t203 = (-t130 * t201 + t134 * t202) / t127 ^ 2;
	t143 = -t171 * t181 + t180 * t193;
	t200 = t134 * t143;
	t188 = t154 * t162 * t170;
	t142 = -t169 * t181 - t180 * t192;
	t179 = -t153 * t157 + t185;
	t155 = t153 * t154;
	t149 = t182 * qJD(2);
	t125 = 0.1e1 / t127;
	t117 = 0.1e1 / t119;
	t116 = t179 * t199;
	t113 = (t136 * t157 + t137 * t194) * t164 - t184 * t116;
	t111 = t179 / t140 ^ 2 * (-t144 * t155 * t186 + t150 * t198) * t206 + (-t149 * t153 * t164 + (-t150 * t188 + (-t157 * t188 + (t155 * t165 * t170 ^ 2 * t206 - t154 * t196) * t146) * qJD(2)) * t165) * t138;
	t1 = [0, t111, 0, 0, 0, 0; 0, (-(t123 * t116 * t115 + t184 * t111) * t121 * t117 + 0.2e1 * (t117 * t204 + t121 * t205) * t113) * t147 + (0.2e1 * t180 * t120 * t205 + (-t151 * t120 + (t180 * t112 + t113 * t152 + (-(t115 * t157 - t116 * t150 + t172 * t189) * t137 - (-t149 + (qJD(2) * t116 - t115) * t194) * t136) * t147) * t121) * t117) * t164, 0, 0, 0, 0; 0, 0.2e1 * (-t131 * t142 + t132 * t200) * t203 + ((t143 * qJD(3) - t151 * t192 + t152 * t169) * t131 + 0.2e1 * t200 * t201 + (-t142 * t129 - (-t142 * qJD(3) + t151 * t193 + t152 * t171) * t134 - t143 * t128) * t132) * t125, -0.2e1 * t203 + 0.2e1 * (t125 * t202 + (-t125 * t201 - t132 * t203) * t134) * t134, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:27
	% DurationCPUTime: 0.79s
	% Computational Cost: add. (1687->72), mult. (5242->175), div. (282->12), fcn. (6739->15), ass. (0->87)
	t200 = sin(pkin(7));
	t197 = t200 ^ 2;
	t246 = 0.2e1 * t197;
	t198 = sin(pkin(13));
	t202 = cos(pkin(13));
	t206 = sin(qJ(3));
	t208 = cos(qJ(3));
	t192 = t206 * t198 - t208 * t202;
	t204 = cos(pkin(7));
	t177 = t192 * t204;
	t199 = sin(pkin(12));
	t203 = cos(pkin(12));
	t209 = cos(qJ(2));
	t205 = cos(pkin(6));
	t207 = sin(qJ(2));
	t234 = t205 * t207;
	t218 = t199 * t234 - t203 * t209;
	t233 = t205 * t209;
	t219 = t199 * t233 + t203 * t207;
	t221 = t208 * t198 + t206 * t202;
	t201 = sin(pkin(6));
	t225 = t199 * t200 * t201;
	t158 = t177 * t219 - t192 * t225 + t218 * t221;
	t153 = t158 ^ 2;
	t178 = t221 * t204;
	t216 = -t178 * t219 + t192 * t218 + t221 * t225;
	t155 = 0.1e1 / t216 ^ 2;
	t245 = t153 * t155;
	t220 = -t199 * t207 + t203 * t233;
	t236 = t201 * t204;
	t171 = t220 * t200 + t203 * t236;
	t237 = t200 * t209;
	t186 = -t201 * t237 + t205 * t204;
	t168 = atan2(t171, t186);
	t163 = sin(t168);
	t164 = cos(t168);
	t151 = t163 * t171 + t164 * t186;
	t148 = 0.1e1 / t151;
	t154 = 0.1e1 / t216;
	t183 = 0.1e1 / t186;
	t149 = 0.1e1 / t151 ^ 2;
	t184 = 0.1e1 / t186 ^ 2;
	t172 = t199 * t236 + t200 * t219;
	t170 = t172 ^ 2;
	t145 = t170 * t149 + 0.1e1;
	t182 = t218 * qJD(2);
	t187 = -t199 * t209 - t203 * t234;
	t180 = t187 * qJD(2);
	t235 = t201 * t207;
	t239 = t171 * t184;
	t223 = t235 * t239;
	t169 = t171 ^ 2;
	t167 = t169 * t184 + 0.1e1;
	t165 = 0.1e1 / t167;
	t240 = t165 * t200;
	t138 = (-qJD(2) * t223 + t180 * t183) * t240;
	t222 = -t163 * t186 + t164 * t171;
	t230 = qJD(2) * t201;
	t224 = t207 * t230;
	t136 = (t163 * t180 + t164 * t224) * t200 + t222 * t138;
	t243 = t136 * t148 * t149;
	t244 = (-t172 * t149 * t182 * t200 - t170 * t243) / t145 ^ 2;
	t174 = qJD(3) * t177;
	t181 = t219 * qJD(2);
	t190 = t192 * qJD(3);
	t191 = t221 * qJD(3);
	t147 = t174 * t219 + t182 * t178 + t181 * t192 - t190 * t225 + t191 * t218;
	t156 = t154 * t155;
	t242 = t156 * t147;
	t162 = t178 * t218 + t192 * t219;
	t241 = t158 * t162;
	t142 = 0.1e1 + t245;
	t175 = t204 * t191;
	t146 = t175 * t219 - t182 * t177 + t181 * t221 - t190 * t218 - t191 * t225;
	t227 = t158 * t155 * t146;
	t228 = 0.2e1 * (-t153 * t242 + t227) / t142 ^ 2;
	t226 = t184 * t197 * t207;
	t217 = -t183 * t187 + t223;
	t185 = t183 * t184;
	t179 = t220 * qJD(2);
	t161 = t177 * t218 - t219 * t221;
	t143 = 0.1e1 / t145;
	t140 = 0.1e1 / t142;
	t139 = t217 * t240;
	t137 = (t163 * t187 + t164 * t235) * t200 - t222 * t139;
	t134 = t217 / t167 ^ 2 * (-t169 * t185 * t224 + t180 * t239) * t246 + (-t179 * t183 * t200 + (-t180 * t226 + (-t187 * t226 + (t185 * t201 * t207 ^ 2 * t246 - t184 * t237) * t171) * qJD(2)) * t201) * t165;
	t1 = [0, t134, 0, 0, 0, 0; 0, (-(t151 * t139 * t138 + t222 * t134) * t149 * t143 + 0.2e1 * (t143 * t243 + t149 * t244) * t137) * t172 + (0.2e1 * t218 * t148 * t244 + (-t181 * t148 + (t218 * t136 + t137 * t182 + (-(t138 * t187 - t139 * t180 + t209 * t230) * t164 - (-t179 + (qJD(2) * t139 - t138) * t235) * t163) * t172) * t149) * t143) * t200, 0, 0, 0, 0; 0, (-t154 * t161 - t155 * t241) * t228 + ((t175 * t218 + t181 * t177 + t182 * t221 + t190 * t219) * t154 - 0.2e1 * t241 * t242 + (-t161 * t147 + (-t174 * t218 + t181 * t178 - t182 * t192 + t191 * t219) * t158 + t162 * t146) * t155) * t140, (-t154 * t216 - t245) * t228 + (0.2e1 * t227 + (-0.2e1 * t153 * t156 - t155 * t216 + t154) * t147) * t140, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:27
	% EndTime: 2019-10-09 22:29:29
	% DurationCPUTime: 2.60s
	% Computational Cost: add. (9541->144), mult. (28576->272), div. (538->12), fcn. (36806->17), ass. (0->121)
	t319 = sin(pkin(13));
	t326 = sin(qJ(3));
	t385 = cos(pkin(13));
	t387 = cos(qJ(3));
	t312 = -t387 * t319 - t326 * t385;
	t389 = qJD(3) * t312;
	t323 = cos(pkin(7));
	t342 = -t326 * t319 + t387 * t385;
	t299 = t342 * t323;
	t327 = sin(qJ(2));
	t329 = cos(qJ(2));
	t388 = t329 * t299 + t327 * t312;
	t324 = cos(pkin(6));
	t386 = cos(pkin(12));
	t356 = t386 * t329;
	t352 = t324 * t356;
	t320 = sin(pkin(12));
	t374 = t320 * t327;
	t305 = t352 - t374;
	t322 = sin(pkin(6));
	t321 = sin(pkin(7));
	t340 = t342 * t321;
	t338 = t322 * t340;
	t357 = t386 * t327;
	t373 = t320 * t329;
	t343 = -t324 * t357 - t373;
	t276 = t305 * t299 - t312 * t343 - t338 * t386;
	t285 = -t388 * t322 - t324 * t340;
	t259 = atan2(t276, t285);
	t254 = sin(t259);
	t255 = cos(t259);
	t244 = t254 * t276 + t255 * t285;
	t241 = 0.1e1 / t244;
	t344 = t324 * t373 + t357;
	t375 = t320 * t322;
	t291 = t321 * t344 + t323 * t375;
	t325 = sin(qJ(5));
	t328 = cos(qJ(5));
	t298 = t312 * t321;
	t300 = t312 * t323;
	t345 = t324 * t374 - t356;
	t341 = -t298 * t375 + t300 * t344 - t342 * t345;
	t267 = t291 * t325 + t328 * t341;
	t263 = 0.1e1 / t267;
	t281 = 0.1e1 / t285;
	t242 = 0.1e1 / t244 ^ 2;
	t264 = 0.1e1 / t267 ^ 2;
	t282 = 0.1e1 / t285 ^ 2;
	t297 = t323 * t389;
	t366 = qJD(2) * t327;
	t301 = -qJD(2) * t352 + t320 * t366;
	t302 = t343 * qJD(2);
	t309 = t342 * qJD(3);
	t337 = t321 * t389;
	t336 = t322 * t337;
	t251 = t305 * t297 + t302 * t299 - t301 * t312 + t309 * t343 - t336 * t386;
	t273 = t276 ^ 2;
	t258 = t273 * t282 + 0.1e1;
	t256 = 0.1e1 / t258;
	t269 = -t324 * t337 + ((-qJD(2) * t312 - t297) * t329 + t299 * t366 + t327 * t309) * t322;
	t376 = t276 * t282;
	t234 = (t251 * t281 - t269 * t376) * t256;
	t350 = -t254 * t285 + t255 * t276;
	t230 = t234 * t350 + t254 * t251 + t255 * t269;
	t384 = t230 * t241 * t242;
	t295 = t321 * t309;
	t296 = qJD(3) * t299;
	t303 = t344 * qJD(2);
	t304 = t345 * qJD(2);
	t253 = t295 * t375 - t296 * t344 - t304 * t300 - t303 * t342 - t345 * t389;
	t371 = t321 * t328;
	t245 = qJD(5) * t267 + t253 * t325 + t304 * t371;
	t266 = -t291 * t328 + t325 * t341;
	t262 = t266 ^ 2;
	t249 = t262 * t264 + 0.1e1;
	t378 = t264 * t266;
	t364 = qJD(5) * t266;
	t372 = t321 * t325;
	t246 = t253 * t328 - t304 * t372 - t364;
	t381 = t246 * t263 * t264;
	t383 = (t245 * t378 - t262 * t381) / t249 ^ 2;
	t278 = -t299 * t344 - t312 * t345 + t320 * t338;
	t382 = t242 * t278;
	t380 = t254 * t278;
	t379 = t255 * t278;
	t377 = t269 * t281 * t282;
	t274 = t278 ^ 2;
	t240 = t274 * t242 + 0.1e1;
	t252 = -t297 * t344 + t304 * t299 - t303 * t312 + t309 * t345 + t320 * t336;
	t363 = 0.2e1 * (t252 * t382 - t274 * t384) / t240 ^ 2;
	t362 = 0.2e1 * t383;
	t361 = 0.2e1 * (t251 * t376 - t273 * t377) / t258 ^ 2;
	t359 = t322 * t386;
	t355 = 0.2e1 * t266 * t381;
	t354 = 0.2e1 * t276 * t377;
	t353 = -0.2e1 * t278 * t384;
	t349 = -t325 * t263 + t328 * t378;
	t275 = -t298 * t359 + t305 * t300 + t342 * t343;
	t284 = -t324 * t298 + (-t300 * t329 + t327 * t342) * t322;
	t348 = -t275 * t281 + t284 * t376;
	t286 = t299 * t343 + t305 * t312;
	t289 = (t299 * t327 - t312 * t329) * t322;
	t347 = -t281 * t286 + t289 * t376;
	t288 = -t300 * t345 - t342 * t344;
	t346 = -t288 * t325 - t345 * t371;
	t272 = t288 * t328 - t345 * t372;
	t287 = -t299 * t345 + t312 * t344;
	t270 = (t388 * qJD(2) + t297 * t327 + t309 * t329) * t322;
	t268 = t324 * t295 + (t296 * t329 + t389 * t327 + (t300 * t327 + t329 * t342) * qJD(2)) * t322;
	t261 = t296 * t345 - t303 * t300 + t304 * t342 - t344 * t389;
	t260 = t297 * t343 + t301 * t299 + t302 * t312 - t305 * t309;
	t250 = t295 * t359 - t305 * t296 + t302 * t300 + t301 * t342 + t343 * t389;
	t247 = 0.1e1 / t249;
	t238 = 0.1e1 / t240;
	t236 = t347 * t256;
	t235 = t348 * t256;
	t232 = -t236 * t350 + t254 * t286 + t255 * t289;
	t231 = -t235 * t350 + t254 * t275 + t255 * t284;
	t229 = t347 * t361 + (t289 * t354 + t260 * t281 + (-t251 * t289 - t269 * t286 - t270 * t276) * t282) * t256;
	t227 = t348 * t361 + (t284 * t354 + t250 * t281 + (-t251 * t284 - t268 * t276 - t269 * t275) * t282) * t256;
	t1 = [0, t229, t227, 0, 0, 0; 0, (-t232 * t382 - t241 * t287) * t363 + ((-t297 * t345 - t303 * t299 - t304 * t312 - t309 * t344) * t241 + t232 * t353 + (-t287 * t230 + t232 * t252 + (t229 * t276 - t236 * t251 + t270 + (t236 * t285 + t286) * t234) * t379 + (-t229 * t285 + t236 * t269 + t260 + (t236 * t276 - t289) * t234) * t380) * t242) * t238, (-t231 * t382 - t241 * t341) * t363 + (t253 * t241 + t231 * t353 + (-t341 * t230 + t231 * t252 + (t227 * t276 - t235 * t251 + t268 + (t235 * t285 + t275) * t234) * t379 + (-t227 * t285 + t235 * t269 + t250 + (t235 * t276 - t284) * t234) * t380) * t242) * t238, 0, 0, 0; 0, (t263 * t346 + t272 * t378) * t362 + ((qJD(5) * t272 + t261 * t325 + t303 * t371) * t263 + t272 * t355 + (t346 * t246 - (qJD(5) * t346 + t261 * t328 - t303 * t372) * t266 - t272 * t245) * t264) * t247, t349 * t278 * t362 + (-t349 * t252 + ((qJD(5) * t263 + t355) * t328 + (-t245 * t328 + (-t246 + t364) * t325) * t264) * t278) * t247, 0, -0.2e1 * t383 + 0.2e1 * (t245 * t264 * t247 + (-t247 * t381 - t264 * t383) * t266) * t266, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:28
	% EndTime: 2019-10-09 22:29:33
	% DurationCPUTime: 5.38s
	% Computational Cost: add. (22420->219), mult. (65370->407), div. (816->12), fcn. (84818->19), ass. (0->159)
	t429 = sin(pkin(7));
	t428 = sin(pkin(13));
	t433 = sin(qJ(3));
	t504 = cos(pkin(13));
	t508 = cos(qJ(3));
	t453 = -t433 * t428 + t508 * t504;
	t407 = t453 * t429;
	t403 = qJD(3) * t407;
	t506 = cos(pkin(7));
	t464 = t506 * t504;
	t477 = t428 * t506;
	t409 = -t433 * t477 + t508 * t464;
	t405 = t409 * qJD(3);
	t410 = t433 * t464 + t508 * t477;
	t434 = sin(qJ(2));
	t437 = cos(qJ(2));
	t505 = cos(pkin(12));
	t507 = cos(pkin(6));
	t466 = t507 * t505;
	t503 = sin(pkin(12));
	t449 = t503 * t434 - t437 * t466;
	t411 = t449 * qJD(2);
	t448 = -t434 * t466 - t503 * t437;
	t412 = t448 * qJD(2);
	t422 = -t508 * t428 - t433 * t504;
	t420 = t422 * qJD(3);
	t430 = sin(pkin(6));
	t476 = t430 * t505;
	t363 = -t403 * t476 - t449 * t405 + t412 * t410 - t411 * t453 - t420 * t448;
	t408 = t422 * t429;
	t388 = t408 * t476 - t449 * t410 - t448 * t453;
	t432 = sin(qJ(5));
	t436 = cos(qJ(5));
	t444 = t449 * t429 - t506 * t476;
	t375 = t388 * t436 + t444 * t432;
	t489 = t429 * t436;
	t339 = t375 * qJD(5) + t363 * t432 + t412 * t489;
	t461 = -t410 * t434 + t437 * t453;
	t379 = t507 * t403 + (t461 * qJD(2) + t405 * t437 + t420 * t434) * t430;
	t462 = t410 * t437 + t434 * t453;
	t395 = -t507 * t408 + t462 * t430;
	t488 = t429 * t437;
	t415 = -t430 * t488 + t507 * t506;
	t393 = t395 * t436 + t415 * t432;
	t478 = t429 * t430 * t434;
	t470 = qJD(2) * t478;
	t345 = t393 * qJD(5) + t379 * t432 - t436 * t470;
	t373 = t388 * t432 - t444 * t436;
	t371 = t373 ^ 2;
	t392 = t395 * t432 - t415 * t436;
	t385 = 0.1e1 / t392 ^ 2;
	t357 = t371 * t385 + 0.1e1;
	t355 = 0.1e1 / t357;
	t384 = 0.1e1 / t392;
	t493 = t373 * t385;
	t322 = (-t339 * t384 + t345 * t493) * t355;
	t358 = atan2(-t373, t392);
	t353 = sin(t358);
	t354 = cos(t358);
	t463 = -t353 * t392 - t354 * t373;
	t317 = t463 * t322 - t353 * t339 + t354 * t345;
	t335 = -t353 * t373 + t354 * t392;
	t332 = 0.1e1 / t335;
	t333 = 0.1e1 / t335 ^ 2;
	t511 = t317 * t332 * t333;
	t465 = t507 * t503;
	t446 = t505 * t434 + t437 * t465;
	t447 = t434 * t465 - t505 * t437;
	t475 = t430 * t503;
	t450 = -t408 * t475 - t410 * t446 - t447 * t453;
	t452 = t429 * t446 + t506 * t475;
	t376 = t432 * t450 - t452 * t436;
	t473 = 0.2e1 * t376 * t511;
	t387 = -t407 * t476 - t449 * t409 - t422 * t448;
	t394 = t507 * t407 + (t409 * t437 + t422 * t434) * t430;
	t456 = -t384 * t387 + t394 * t493;
	t510 = t432 * t456;
	t499 = t345 * t384 * t385;
	t509 = -0.2e1 * (t339 * t493 - t371 * t499) / t357 ^ 2;
	t377 = t452 * t432 + t436 * t450;
	t390 = t407 * t475 - t409 * t446 - t422 * t447;
	t431 = sin(qJ(6));
	t435 = cos(qJ(6));
	t352 = t377 * t435 - t390 * t431;
	t348 = 0.1e1 / t352;
	t349 = 0.1e1 / t352 ^ 2;
	t414 = t447 * qJD(2);
	t413 = t446 * qJD(2);
	t445 = t403 * t475 - t405 * t446 + t414 * t410 - t413 * t453 - t420 * t447;
	t490 = t429 * t432;
	t342 = -t376 * qJD(5) - t414 * t490 + t436 * t445;
	t404 = t429 * t420;
	t406 = t410 * qJD(3);
	t419 = t453 * qJD(3);
	t365 = t404 * t475 + t406 * t446 + t414 * t409 - t413 * t422 + t419 * t447;
	t351 = t377 * t431 + t390 * t435;
	t484 = qJD(6) * t351;
	t331 = t342 * t435 - t365 * t431 - t484;
	t502 = t331 * t348 * t349;
	t501 = t333 * t376;
	t330 = t352 * qJD(6) + t342 * t431 + t365 * t435;
	t347 = t351 ^ 2;
	t338 = t347 * t349 + 0.1e1;
	t497 = t349 * t351;
	t500 = 0.1e1 / t338 ^ 2 * (t330 * t497 - t347 * t502);
	t498 = t348 * t431;
	t496 = t351 * t435;
	t495 = t353 * t376;
	t494 = t354 * t376;
	t492 = t390 * t432;
	t491 = t390 * t436;
	t485 = qJD(5) * t436;
	t372 = t376 ^ 2;
	t329 = t333 * t372 + 0.1e1;
	t341 = t377 * qJD(5) + t414 * t489 + t432 * t445;
	t483 = 0.2e1 * (t341 * t501 - t372 * t511) / t329 ^ 2;
	t481 = -0.2e1 * t500;
	t480 = 0.2e1 * t500;
	t479 = t351 * t502;
	t472 = 0.2e1 * t479;
	t471 = -0.2e1 * t373 * t499;
	t467 = qJD(6) * t491 - t445;
	t398 = t410 * t447 - t446 * t453;
	t383 = t398 * t436 - t447 * t490;
	t397 = -t409 * t447 + t422 * t446;
	t368 = t383 * t435 + t397 * t431;
	t367 = t383 * t431 - t397 * t435;
	t459 = t349 * t496 - t498;
	t458 = -t375 * t384 + t393 * t493;
	t396 = t410 * t448 - t449 * t453;
	t381 = t396 * t432 + t448 * t489;
	t400 = t461 * t430;
	t399 = t400 * t432 - t436 * t478;
	t457 = -t381 * t384 + t399 * t493;
	t455 = -t398 * t432 - t447 * t489;
	t451 = -qJD(5) * t492 + qJD(6) * t450 + t365 * t436;
	t378 = t507 * t404 + (-t406 * t437 - t419 * t434 + (-t409 * t434 + t422 * t437) * qJD(2)) * t430;
	t370 = t405 * t447 + t410 * t413 + t414 * t453 - t420 * t446;
	t369 = t406 * t447 - t409 * t413 - t414 * t422 - t419 * t446;
	t362 = -t404 * t476 + t449 * t406 + t412 * t409 - t411 * t422 + t419 * t448;
	t361 = t400 * t485 + ((t420 * t437 + (qJD(5) * t429 - t405) * t434) * t432 + (-t462 * t432 - t436 * t488) * qJD(2)) * t430;
	t360 = t431 * t450 + t435 * t491;
	t359 = t431 * t491 - t435 * t450;
	t346 = -t392 * qJD(5) + t379 * t436 + t432 * t470;
	t344 = t455 * qJD(5) + t370 * t436 - t413 * t490;
	t343 = (t405 * t448 + t411 * t410 + t412 * t453 - t449 * t420) * t432 + t411 * t489 + (t396 * t436 - t448 * t490) * qJD(5);
	t340 = -t373 * qJD(5) + t363 * t436 - t412 * t490;
	t336 = 0.1e1 / t338;
	t327 = 0.1e1 / t329;
	t326 = t457 * t355;
	t325 = t355 * t510;
	t324 = t458 * t355;
	t320 = t463 * t326 - t353 * t381 + t354 * t399;
	t319 = (-t353 * t387 + t354 * t394) * t432 + t463 * t325;
	t318 = t463 * t324 - t353 * t375 + t354 * t393;
	t316 = t457 * t509 + (t399 * t471 - t343 * t384 + (t339 * t399 + t345 * t381 + t361 * t373) * t385) * t355;
	t314 = t458 * t509 + (t393 * t471 - t340 * t384 + (t339 * t393 + t345 * t375 + t346 * t373) * t385) * t355;
	t313 = t509 * t510 + (t456 * t485 + (t394 * t471 - t362 * t384 + (t339 * t394 + t345 * t387 + t373 * t378) * t385) * t432) * t355;
	t1 = [0, t316, t313, 0, t314, 0; 0, (t320 * t501 + t332 * t455) * t483 + ((t383 * qJD(5) + t370 * t432 + t413 * t489) * t332 + t320 * t473 + (t455 * t317 - t320 * t341 - (-t316 * t373 - t326 * t339 + t361 + (-t326 * t392 - t381) * t322) * t494 - (-t316 * t392 - t326 * t345 - t343 + (t326 * t373 - t399) * t322) * t495) * t333) * t327, (t319 * t501 - t332 * t492) * t483 + ((t365 * t432 + t390 * t485) * t332 + t319 * t473 + (-t319 * t341 - t492 * t317 - (t394 * t485 - t313 * t373 - t325 * t339 + t378 * t432 + (-t325 * t392 - t387 * t432) * t322) * t494 - (-t387 * t485 - t313 * t392 - t325 * t345 - t362 * t432 + (t325 * t373 - t394 * t432) * t322) * t495) * t333) * t327, 0, (t318 * t501 - t332 * t377) * t483 + (t318 * t473 + t342 * t332 + (-t377 * t317 - t318 * t341 - (-t314 * t373 - t324 * t339 + t346 + (-t324 * t392 - t375) * t322) * t494 - (-t314 * t392 - t324 * t345 - t340 + (t324 * t373 - t393) * t322) * t495) * t333) * t327, 0; 0, (-t348 * t367 + t368 * t497) * t480 + ((t368 * qJD(6) + t344 * t431 - t369 * t435) * t348 + t368 * t472 + (-t367 * t331 - (-t367 * qJD(6) + t344 * t435 + t369 * t431) * t351 - t368 * t330) * t349) * t336, (-t348 * t359 + t360 * t497) * t480 + (t360 * t472 + t467 * t348 * t435 + t451 * t498 + (t467 * t351 * t431 - t360 * t330 - t359 * t331 - t451 * t496) * t349) * t336, 0, t459 * t376 * t481 + (t459 * t341 + ((-qJD(6) * t348 - 0.2e1 * t479) * t435 + (t330 * t435 + (t331 - t484) * t431) * t349) * t376) * t336, t481 + 0.2e1 * (t330 * t349 * t336 + (-t336 * t502 - t349 * t500) * t351) * t351;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end