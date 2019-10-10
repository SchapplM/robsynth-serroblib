% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR1
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
%   Wie in S6PRPRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(11));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(11)) * cos(pkin(6));
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
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (153->13), mult. (451->28), div. (25->4), fcn. (552->7), ass. (0->20)
	t63 = sin(pkin(12));
	t64 = cos(pkin(12));
	t65 = cos(pkin(11));
	t66 = sin(qJ(2));
	t67 = cos(qJ(2));
	t78 = cos(pkin(6));
	t74 = t67 * t78;
	t75 = t66 * t78;
	t77 = sin(pkin(11));
	t54 = (t63 * t75 - t64 * t74) * t77 + t65 * (-t67 * t63 - t66 * t64);
	t48 = t54 * qJD(2);
	t72 = (-t63 * t74 - t64 * t75) * t77 - t65 * (t66 * t63 - t67 * t64);
	t51 = 0.1e1 / t72 ^ 2;
	t84 = t51 * t54 ^ 2;
	t50 = 0.1e1 / t72;
	t83 = t50 * t84;
	t82 = t51 * t72;
	t76 = t82 * t48;
	t46 = 0.1e1 + t84;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0.2e1 * (-t50 * t72 - t84) / t46 ^ 2 * (-t48 * t83 - t76) + (-0.2e1 * t76 + (t50 - t82 - 0.2e1 * t83) * t48) / t46, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:49
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (1747->58), mult. (5333->135), div. (281->12), fcn. (6885->13), ass. (0->71)
	t185 = sin(qJ(4));
	t187 = cos(qJ(4));
	t184 = cos(pkin(6));
	t179 = sin(pkin(12));
	t182 = cos(pkin(12));
	t186 = sin(qJ(2));
	t188 = cos(qJ(2));
	t200 = t188 * t179 + t186 * t182;
	t171 = t200 * t184;
	t174 = t186 * t179 - t188 * t182;
	t180 = sin(pkin(11));
	t183 = cos(pkin(11));
	t201 = -t180 * t171 - t183 * t174;
	t181 = sin(pkin(6));
	t205 = t180 * t181;
	t197 = -t185 * t201 + t187 * t205;
	t214 = t197 * qJD(4);
	t172 = t174 * qJD(2);
	t196 = t174 * t184;
	t156 = -t180 * t200 - t183 * t196;
	t169 = t174 * t181;
	t142 = atan2(t156, t169);
	t137 = sin(t142);
	t138 = cos(t142);
	t131 = t137 * t156 + t138 * t169;
	t128 = 0.1e1 / t131;
	t148 = t185 * t205 + t187 * t201;
	t144 = 0.1e1 / t148;
	t166 = 0.1e1 / t169;
	t129 = 0.1e1 / t131 ^ 2;
	t145 = 0.1e1 / t148 ^ 2;
	t167 = 0.1e1 / t169 ^ 2;
	t153 = t156 ^ 2;
	t141 = t153 * t167 + 0.1e1;
	t139 = 0.1e1 / t141;
	t195 = qJD(2) * t171;
	t149 = t180 * t172 - t183 * t195;
	t170 = t200 * t181;
	t163 = qJD(2) * t170;
	t208 = t156 * t167;
	t122 = (t149 * t166 - t163 * t208) * t139;
	t202 = -t137 * t169 + t138 * t156;
	t119 = t202 * t122 + t137 * t149 + t138 * t163;
	t213 = t119 * t128 * t129;
	t143 = t197 ^ 2;
	t134 = t143 * t145 + 0.1e1;
	t165 = t184 * t172;
	t173 = t200 * qJD(2);
	t152 = t180 * t165 - t183 * t173;
	t135 = t148 * qJD(4) + t152 * t185;
	t209 = t145 * t197;
	t136 = t152 * t187 + t214;
	t210 = t136 * t144 * t145;
	t212 = (-t135 * t209 - t143 * t210) / t134 ^ 2;
	t158 = t180 * t196 - t183 * t200;
	t211 = t129 * t158;
	t207 = t156 * t170;
	t206 = t163 * t166 * t167;
	t199 = -t144 * t185 - t187 * t209;
	t155 = -t183 * t171 + t180 * t174;
	t198 = -t155 * t166 + t167 * t207;
	t164 = t181 * t172;
	t154 = t158 ^ 2;
	t151 = t183 * t172 + t180 * t195;
	t150 = t183 * t165 + t180 * t173;
	t132 = 0.1e1 / t134;
	t126 = t154 * t129 + 0.1e1;
	t123 = t198 * t139;
	t120 = -t202 * t123 + t137 * t155 + t138 * t170;
	t118 = 0.2e1 * t198 / t141 ^ 2 * (t149 * t208 - t153 * t206) + (0.2e1 * t206 * t207 + t150 * t166 + (-t149 * t170 - t155 * t163 + t156 * t164) * t167) * t139;
	t1 = [0, t118, 0, 0, 0, 0; 0, 0.2e1 * (-t120 * t211 - t128 * t201) / t126 ^ 2 * (t151 * t211 - t154 * t213) + (t152 * t128 + (-t119 * t201 + t120 * t151) * t129 + (-0.2e1 * t120 * t213 + ((t118 * t156 - t123 * t149 - t164 + (t123 * t169 + t155) * t122) * t138 + (-t118 * t169 + t123 * t163 + t150 + (t123 * t156 - t170) * t122) * t137) * t129) * t158) / t126, 0, 0, 0, 0; 0, 0.2e1 * t199 * t158 * t212 + (-t199 * t151 + ((qJD(4) * t144 - 0.2e1 * t197 * t210) * t187 + (-t135 * t187 + (-t136 - t214) * t185) * t145) * t158) * t132, 0, -0.2e1 * t212 - 0.2e1 * (t132 * t135 * t145 - (-t132 * t210 - t145 * t212) * t197) * t197, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:49
	% DurationCPUTime: 0.71s
	% Computational Cost: add. (2232->61), mult. (5863->138), div. (299->12), fcn. (7536->13), ass. (0->75)
	t217 = cos(pkin(6));
	t212 = sin(pkin(12));
	t215 = cos(pkin(12));
	t218 = sin(qJ(2));
	t219 = cos(qJ(2));
	t230 = t219 * t212 + t218 * t215;
	t200 = t230 * t217;
	t203 = t218 * t212 - t219 * t215;
	t201 = t203 * qJD(2);
	t213 = sin(pkin(11));
	t216 = cos(pkin(11));
	t227 = t203 * t217;
	t185 = -t213 * t230 - t216 * t227;
	t214 = sin(pkin(6));
	t198 = t203 * t214;
	t171 = atan2(t185, t198);
	t166 = sin(t171);
	t167 = cos(t171);
	t160 = t166 * t185 + t167 * t198;
	t157 = 0.1e1 / t160;
	t211 = qJ(4) + qJ(5);
	t208 = sin(t211);
	t209 = cos(t211);
	t231 = -t213 * t200 - t216 * t203;
	t236 = t213 * t214;
	t177 = t208 * t236 + t209 * t231;
	t173 = 0.1e1 / t177;
	t195 = 0.1e1 / t198;
	t158 = 0.1e1 / t160 ^ 2;
	t174 = 0.1e1 / t177 ^ 2;
	t196 = 0.1e1 / t198 ^ 2;
	t182 = t185 ^ 2;
	t170 = t182 * t196 + 0.1e1;
	t168 = 0.1e1 / t170;
	t226 = qJD(2) * t200;
	t178 = t213 * t201 - t216 * t226;
	t199 = t230 * t214;
	t192 = qJD(2) * t199;
	t240 = t185 * t196;
	t151 = (t178 * t195 - t192 * t240) * t168;
	t232 = -t166 * t198 + t167 * t185;
	t148 = t232 * t151 + t166 * t178 + t167 * t192;
	t245 = t148 * t157 * t158;
	t176 = t208 * t231 - t209 * t236;
	t172 = t176 ^ 2;
	t163 = t172 * t174 + 0.1e1;
	t194 = t217 * t201;
	t202 = t230 * qJD(2);
	t181 = t213 * t194 - t216 * t202;
	t210 = qJD(4) + qJD(5);
	t233 = t210 * t236 + t181;
	t238 = t231 * t210;
	t164 = t233 * t208 + t209 * t238;
	t241 = t174 * t176;
	t165 = -t208 * t238 + t233 * t209;
	t242 = t165 * t173 * t174;
	t244 = (t164 * t241 - t172 * t242) / t163 ^ 2;
	t187 = t213 * t227 - t216 * t230;
	t243 = t158 * t187;
	t239 = t185 * t199;
	t237 = t192 * t195 * t196;
	t229 = -t173 * t208 + t209 * t241;
	t184 = -t216 * t200 + t213 * t203;
	t228 = -t184 * t195 + t196 * t239;
	t193 = t214 * t201;
	t183 = t187 ^ 2;
	t180 = t216 * t201 + t213 * t226;
	t179 = t216 * t194 + t213 * t202;
	t161 = 0.1e1 / t163;
	t155 = t183 * t158 + 0.1e1;
	t152 = t228 * t168;
	t149 = -t232 * t152 + t166 * t184 + t167 * t199;
	t147 = 0.2e1 * t228 / t170 ^ 2 * (t178 * t240 - t182 * t237) + (0.2e1 * t237 * t239 + t179 * t195 + (-t178 * t199 - t184 * t192 + t185 * t193) * t196) * t168;
	t145 = -0.2e1 * t244 + 0.2e1 * (t161 * t164 * t174 + (-t161 * t242 - t174 * t244) * t176) * t176;
	t1 = [0, t147, 0, 0, 0, 0; 0, 0.2e1 * (-t149 * t243 - t157 * t231) / t155 ^ 2 * (t180 * t243 - t183 * t245) + (t181 * t157 + (-t148 * t231 + t149 * t180) * t158 + (-0.2e1 * t149 * t245 + ((t147 * t185 - t152 * t178 - t193 + (t152 * t198 + t184) * t151) * t167 + (-t147 * t198 + t152 * t192 + t179 + (t152 * t185 - t199) * t151) * t166) * t158) * t187) / t155, 0, 0, 0, 0; 0, 0.2e1 * t229 * t187 * t244 + (-t229 * t180 + ((t173 * t210 + 0.2e1 * t176 * t242) * t209 + (-t164 * t209 + (t176 * t210 - t165) * t208) * t174) * t187) * t161, 0, t145, t145, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:49
	% EndTime: 2019-10-09 21:53:51
	% DurationCPUTime: 2.35s
	% Computational Cost: add. (12871->119), mult. (23866->239), div. (822->12), fcn. (31086->15), ass. (0->117)
	t319 = sin(pkin(12));
	t322 = cos(pkin(12));
	t326 = sin(qJ(2));
	t328 = cos(qJ(2));
	t309 = t326 * t319 - t328 * t322;
	t324 = cos(pkin(6));
	t337 = t309 * t324;
	t302 = qJD(2) * t337;
	t341 = t328 * t319 + t326 * t322;
	t307 = t341 * qJD(2);
	t320 = sin(pkin(11));
	t323 = cos(pkin(11));
	t283 = -t323 * t302 - t320 * t307;
	t305 = t341 * t324;
	t289 = t323 * t305 - t320 * t309;
	t318 = qJ(4) + qJ(5);
	t315 = sin(t318);
	t317 = qJD(4) + qJD(5);
	t321 = sin(pkin(6));
	t358 = t321 * t323;
	t349 = t315 * t358;
	t316 = cos(t318);
	t360 = t316 * t317;
	t255 = t283 * t315 + t289 * t360 - t317 * t349;
	t277 = t289 * t315 + t316 * t358;
	t275 = t277 ^ 2;
	t304 = t341 * t321;
	t296 = t304 * t315 - t324 * t316;
	t294 = 0.1e1 / t296 ^ 2;
	t269 = t275 * t294 + 0.1e1;
	t267 = 0.1e1 / t269;
	t303 = t309 * t321;
	t347 = -qJD(2) * t303 + t317 * t324;
	t273 = t304 * t360 + t347 * t315;
	t293 = 0.1e1 / t296;
	t365 = t277 * t294;
	t239 = (-t255 * t293 + t273 * t365) * t267;
	t272 = atan2(-t277, t296);
	t265 = sin(t272);
	t266 = cos(t272);
	t344 = -t265 * t296 - t266 * t277;
	t235 = t344 * t239 - t265 * t255 + t266 * t273;
	t249 = -t265 * t277 + t266 * t296;
	t246 = 0.1e1 / t249;
	t247 = 0.1e1 / t249 ^ 2;
	t379 = t235 * t246 * t247;
	t342 = -t320 * t305 - t323 * t309;
	t359 = t320 * t321;
	t280 = t315 * t342 - t316 * t359;
	t378 = 0.2e1 * t280 * t379;
	t288 = -t320 * t341 - t323 * t337;
	t338 = -t288 * t293 - t303 * t365;
	t377 = t315 * t338;
	t366 = t273 * t293 * t294;
	t376 = -0.2e1 * (t255 * t365 - t275 * t366) / t269 ^ 2;
	t281 = t315 * t359 + t316 * t342;
	t327 = cos(qJ(6));
	t291 = t320 * t337 - t323 * t341;
	t325 = sin(qJ(6));
	t363 = t291 * t325;
	t264 = t281 * t327 - t363;
	t260 = 0.1e1 / t264;
	t261 = 0.1e1 / t264 ^ 2;
	t343 = t320 * t302 - t323 * t307;
	t346 = t317 * t359 + t343;
	t361 = t315 * t317;
	t258 = t346 * t316 - t342 * t361;
	t306 = t309 * qJD(2);
	t336 = t324 * t307;
	t284 = t323 * t306 + t320 * t336;
	t250 = t264 * qJD(6) + t258 * t325 + t284 * t327;
	t362 = t291 * t327;
	t263 = t281 * t325 + t362;
	t259 = t263 ^ 2;
	t254 = t259 * t261 + 0.1e1;
	t370 = t261 * t263;
	t354 = qJD(6) * t263;
	t251 = t258 * t327 - t284 * t325 - t354;
	t373 = t251 * t260 * t261;
	t375 = (t250 * t370 - t259 * t373) / t254 ^ 2;
	t374 = t247 * t280;
	t257 = t346 * t315 + t342 * t360;
	t372 = t257 * t247;
	t371 = t260 * t325;
	t369 = t263 * t327;
	t368 = t265 * t280;
	t367 = t266 * t280;
	t364 = t291 * t315;
	t276 = t280 ^ 2;
	t245 = t276 * t247 + 0.1e1;
	t353 = 0.2e1 * (-t276 * t379 + t280 * t372) / t245 ^ 2;
	t352 = -0.2e1 * t375;
	t350 = t263 * t373;
	t348 = -0.2e1 * t277 * t366;
	t345 = qJD(6) * t291 * t316 - t343;
	t340 = t261 * t369 - t371;
	t279 = t289 * t316 - t349;
	t297 = t304 * t316 + t324 * t315;
	t339 = -t279 * t293 + t297 * t365;
	t335 = qJD(6) * t342 + t284 * t316 - t291 * t361;
	t300 = t321 * t307;
	t282 = t320 * t306 - t323 * t336;
	t274 = -t304 * t361 + t347 * t316;
	t271 = t316 * t362 + t325 * t342;
	t270 = t316 * t363 - t327 * t342;
	t256 = -t289 * t361 + (-t317 * t358 + t283) * t316;
	t252 = 0.1e1 / t254;
	t243 = 0.1e1 / t245;
	t241 = t267 * t377;
	t240 = t339 * t267;
	t237 = (-t265 * t288 - t266 * t303) * t315 + t344 * t241;
	t236 = t344 * t240 - t265 * t279 + t266 * t297;
	t233 = t339 * t376 + (t297 * t348 - t256 * t293 + (t255 * t297 + t273 * t279 + t274 * t277) * t294) * t267;
	t232 = t376 * t377 + (t338 * t360 + (-t303 * t348 - t282 * t293 + (-t255 * t303 + t273 * t288 - t277 * t300) * t294) * t315) * t267;
	t231 = t340 * t280 * t352 + (t340 * t257 + ((-qJD(6) * t260 - 0.2e1 * t350) * t327 + (t250 * t327 + (t251 - t354) * t325) * t261) * t280) * t252;
	t230 = (t236 * t374 - t246 * t281) * t353 + (t236 * t378 + t258 * t246 + (-t281 * t235 - t236 * t257 - (-t233 * t277 - t240 * t255 + t274 + (-t240 * t296 - t279) * t239) * t367 - (-t233 * t296 - t240 * t273 - t256 + (t240 * t277 - t297) * t239) * t368) * t247) * t243;
	t1 = [0, t232, 0, t233, t233, 0; 0, (t237 * t374 - t246 * t364) * t353 + ((t284 * t315 + t291 * t360) * t246 + (-t372 + t378) * t237 + (-t364 * t235 - (-t303 * t360 - t232 * t277 - t241 * t255 - t300 * t315 + (-t241 * t296 - t288 * t315) * t239) * t367 - (-t288 * t360 - t232 * t296 - t241 * t273 - t282 * t315 + (t241 * t277 + t303 * t315) * t239) * t368) * t247) * t243, 0, t230, t230, 0; 0, 0.2e1 * (-t260 * t270 + t271 * t370) * t375 + (0.2e1 * t271 * t350 + t345 * t260 * t327 + t335 * t371 + (t345 * t263 * t325 - t271 * t250 - t270 * t251 - t335 * t369) * t261) * t252, 0, t231, t231, t352 + 0.2e1 * (t250 * t261 * t252 + (-t252 * t373 - t261 * t375) * t263) * t263;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end