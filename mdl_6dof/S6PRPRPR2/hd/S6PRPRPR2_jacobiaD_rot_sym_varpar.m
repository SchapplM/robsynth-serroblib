% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRPR2
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
%   Wie in S6PRPRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(10));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(10)) * cos(pkin(6));
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
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:51
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (153->13), mult. (451->28), div. (25->4), fcn. (552->7), ass. (0->20)
	t63 = sin(pkin(11));
	t64 = cos(pkin(11));
	t65 = cos(pkin(10));
	t66 = sin(qJ(2));
	t67 = cos(qJ(2));
	t78 = cos(pkin(6));
	t74 = t67 * t78;
	t75 = t66 * t78;
	t77 = sin(pkin(10));
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
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:52
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (1747->58), mult. (5333->135), div. (281->12), fcn. (6885->13), ass. (0->71)
	t185 = sin(qJ(4));
	t187 = cos(qJ(4));
	t184 = cos(pkin(6));
	t179 = sin(pkin(11));
	t182 = cos(pkin(11));
	t186 = sin(qJ(2));
	t188 = cos(qJ(2));
	t200 = t188 * t179 + t186 * t182;
	t171 = t200 * t184;
	t174 = t186 * t179 - t188 * t182;
	t180 = sin(pkin(10));
	t183 = cos(pkin(10));
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
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:53
	% DurationCPUTime: 1.72s
	% Computational Cost: add. (5050->106), mult. (14713->222), div. (535->12), fcn. (19181->15), ass. (0->106)
	t244 = sin(pkin(11));
	t248 = cos(pkin(11));
	t252 = sin(qJ(2));
	t254 = cos(qJ(2));
	t236 = t252 * t244 - t254 * t248;
	t250 = cos(pkin(6));
	t262 = t236 * t250;
	t230 = qJD(2) * t262;
	t267 = t254 * t244 + t252 * t248;
	t235 = t267 * qJD(2);
	t245 = sin(pkin(10));
	t249 = cos(pkin(10));
	t212 = -t249 * t230 - t245 * t235;
	t233 = t267 * t250;
	t217 = t249 * t233 - t245 * t236;
	t251 = sin(qJ(4));
	t246 = sin(pkin(6));
	t282 = t246 * t251;
	t272 = t249 * t282;
	t253 = cos(qJ(4));
	t277 = qJD(4) * t253;
	t184 = -qJD(4) * t272 + t212 * t251 + t217 * t277;
	t281 = t246 * t253;
	t206 = t217 * t251 + t249 * t281;
	t204 = t206 ^ 2;
	t232 = t267 * t246;
	t224 = t232 * t251 - t250 * t253;
	t222 = 0.1e1 / t224 ^ 2;
	t200 = t204 * t222 + 0.1e1;
	t198 = 0.1e1 / t200;
	t225 = t232 * t253 + t250 * t251;
	t231 = t236 * t246;
	t229 = qJD(2) * t231;
	t202 = t225 * qJD(4) - t229 * t251;
	t221 = 0.1e1 / t224;
	t285 = t206 * t222;
	t168 = (-t184 * t221 + t202 * t285) * t198;
	t201 = atan2(-t206, t224);
	t196 = sin(t201);
	t197 = cos(t201);
	t270 = -t196 * t224 - t197 * t206;
	t164 = t270 * t168 - t196 * t184 + t197 * t202;
	t178 = -t196 * t206 + t197 * t224;
	t175 = 0.1e1 / t178;
	t176 = 0.1e1 / t178 ^ 2;
	t299 = t164 * t175 * t176;
	t268 = -t245 * t233 - t249 * t236;
	t210 = t245 * t282 + t253 * t268;
	t219 = t245 * t262 - t249 * t267;
	t243 = sin(pkin(12));
	t247 = cos(pkin(12));
	t192 = t210 * t243 + t219 * t247;
	t264 = t245 * t281 - t251 * t268;
	t269 = t245 * t230 - t249 * t235;
	t187 = t264 * qJD(4) + t253 * t269;
	t234 = t236 * qJD(2);
	t261 = t250 * t235;
	t213 = t249 * t234 + t245 * t261;
	t183 = t187 * t247 - t213 * t243;
	t193 = t210 * t247 - t219 * t243;
	t189 = 0.1e1 / t193;
	t190 = 0.1e1 / t193 ^ 2;
	t293 = t183 * t189 * t190;
	t298 = 0.2e1 * t192 * t293;
	t297 = -0.2e1 * t264 * t299;
	t216 = -t245 * t267 - t249 * t262;
	t265 = -t216 * t221 - t231 * t285;
	t296 = t251 * t265;
	t286 = t202 * t221 * t222;
	t295 = -0.2e1 * (t184 * t285 - t204 * t286) / t200 ^ 2;
	t294 = t176 * t264;
	t186 = t210 * qJD(4) + t251 * t269;
	t292 = t186 * t176;
	t291 = t189 * t243;
	t290 = t190 * t192;
	t289 = t192 * t247;
	t288 = t196 * t264;
	t287 = t197 * t264;
	t284 = t219 * t251;
	t283 = t219 * t253;
	t205 = t264 ^ 2;
	t174 = t205 * t176 + 0.1e1;
	t276 = 0.2e1 * (-t205 * t299 - t264 * t292) / t174 ^ 2;
	t188 = t192 ^ 2;
	t181 = t188 * t190 + 0.1e1;
	t182 = t187 * t243 + t213 * t247;
	t275 = 0.2e1 * (t182 * t290 - t188 * t293) / t181 ^ 2;
	t271 = -0.2e1 * t206 * t286;
	t208 = t217 * t253 - t272;
	t266 = -t208 * t221 + t225 * t285;
	t263 = -qJD(4) * t284 + t213 * t253;
	t228 = t246 * t235;
	t211 = t245 * t234 - t249 * t261;
	t203 = -t224 * qJD(4) - t229 * t253;
	t195 = t243 * t268 + t247 * t283;
	t194 = t243 * t283 - t247 * t268;
	t185 = -t206 * qJD(4) + t212 * t253;
	t179 = 0.1e1 / t181;
	t172 = 0.1e1 / t174;
	t170 = t198 * t296;
	t169 = t266 * t198;
	t166 = (-t196 * t216 - t197 * t231) * t251 + t270 * t170;
	t165 = t270 * t169 - t196 * t208 + t197 * t225;
	t163 = t266 * t295 + (t225 * t271 - t185 * t221 + (t184 * t225 + t202 * t208 + t203 * t206) * t222) * t198;
	t161 = t295 * t296 + (t265 * t277 + (-t231 * t271 - t211 * t221 + (-t184 * t231 + t202 * t216 - t206 * t228) * t222) * t251) * t198;
	t1 = [0, t161, 0, t163, 0, 0; 0, (-t166 * t294 - t175 * t284) * t276 + ((t213 * t251 + t219 * t277) * t175 + (-t292 + t297) * t166 + (-t284 * t164 + (-t231 * t277 - t161 * t206 - t170 * t184 - t228 * t251 + (-t170 * t224 - t216 * t251) * t168) * t287 + (-t216 * t277 - t161 * t224 - t170 * t202 - t211 * t251 + (t170 * t206 + t231 * t251) * t168) * t288) * t176) * t172, 0, (-t165 * t294 - t175 * t210) * t276 + (t165 * t297 + t187 * t175 + (-t210 * t164 - t165 * t186 + (-t163 * t206 - t169 * t184 + t203 + (-t169 * t224 - t208) * t168) * t287 + (-t163 * t224 - t169 * t202 - t185 + (t169 * t206 - t225) * t168) * t288) * t176) * t172, 0, 0; 0, (-t189 * t194 + t195 * t290) * t275 + ((t263 * t243 - t247 * t269) * t189 + t195 * t298 + (-t194 * t183 - (t243 * t269 + t263 * t247) * t192 - t195 * t182) * t190) * t179, 0, -(-t190 * t289 + t291) * t264 * t275 + (t264 * t247 * t298 - t186 * t291 + (t186 * t289 - (t182 * t247 + t183 * t243) * t264) * t190) * t179, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:31:51
	% EndTime: 2019-10-09 21:31:53
	% DurationCPUTime: 1.84s
	% Computational Cost: add. (5953->115), mult. (16325->235), div. (559->12), fcn. (21250->15), ass. (0->111)
	t282 = sin(pkin(11));
	t285 = cos(pkin(11));
	t289 = sin(qJ(2));
	t291 = cos(qJ(2));
	t272 = t289 * t282 - t291 * t285;
	t287 = cos(pkin(6));
	t300 = t272 * t287;
	t266 = qJD(2) * t300;
	t305 = t291 * t282 + t289 * t285;
	t271 = t305 * qJD(2);
	t283 = sin(pkin(10));
	t286 = cos(pkin(10));
	t247 = -t286 * t266 - t283 * t271;
	t269 = t305 * t287;
	t253 = t286 * t269 - t283 * t272;
	t288 = sin(qJ(4));
	t284 = sin(pkin(6));
	t322 = t284 * t288;
	t311 = t286 * t322;
	t290 = cos(qJ(4));
	t317 = qJD(4) * t290;
	t225 = -qJD(4) * t311 + t247 * t288 + t253 * t317;
	t321 = t284 * t290;
	t241 = t253 * t288 + t286 * t321;
	t239 = t241 ^ 2;
	t268 = t305 * t284;
	t260 = t268 * t288 - t287 * t290;
	t258 = 0.1e1 / t260 ^ 2;
	t235 = t239 * t258 + 0.1e1;
	t233 = 0.1e1 / t235;
	t261 = t268 * t290 + t287 * t288;
	t267 = t272 * t284;
	t265 = qJD(2) * t267;
	t237 = qJD(4) * t261 - t265 * t288;
	t257 = 0.1e1 / t260;
	t325 = t241 * t258;
	t203 = (-t225 * t257 + t237 * t325) * t233;
	t236 = atan2(-t241, t260);
	t231 = sin(t236);
	t232 = cos(t236);
	t308 = -t231 * t260 - t232 * t241;
	t199 = t203 * t308 - t231 * t225 + t232 * t237;
	t215 = -t231 * t241 + t232 * t260;
	t212 = 0.1e1 / t215;
	t213 = 0.1e1 / t215 ^ 2;
	t339 = t199 * t212 * t213;
	t306 = -t283 * t269 - t286 * t272;
	t301 = t283 * t321 - t288 * t306;
	t338 = -0.2e1 * t301 * t339;
	t252 = -t283 * t305 - t286 * t300;
	t302 = -t252 * t257 - t267 * t325;
	t337 = t288 * t302;
	t326 = t237 * t257 * t258;
	t336 = -0.2e1 * (t225 * t325 - t239 * t326) / t235 ^ 2;
	t245 = t283 * t322 + t290 * t306;
	t255 = t283 * t300 - t286 * t305;
	t281 = pkin(12) + qJ(6);
	t279 = sin(t281);
	t280 = cos(t281);
	t224 = t245 * t280 - t255 * t279;
	t220 = 0.1e1 / t224;
	t221 = 0.1e1 / t224 ^ 2;
	t307 = t283 * t266 - t286 * t271;
	t228 = qJD(4) * t301 + t290 * t307;
	t270 = t272 * qJD(2);
	t299 = t287 * t271;
	t248 = t286 * t270 + t283 * t299;
	t210 = qJD(6) * t224 + t228 * t279 + t248 * t280;
	t223 = t245 * t279 + t255 * t280;
	t219 = t223 ^ 2;
	t218 = t219 * t221 + 0.1e1;
	t331 = t221 * t223;
	t316 = qJD(6) * t223;
	t211 = t228 * t280 - t248 * t279 - t316;
	t334 = t211 * t220 * t221;
	t335 = (t210 * t331 - t219 * t334) / t218 ^ 2;
	t333 = t213 * t301;
	t332 = t220 * t279;
	t330 = t223 * t280;
	t227 = qJD(4) * t245 + t288 * t307;
	t329 = t227 * t213;
	t328 = t231 * t301;
	t327 = t232 * t301;
	t324 = t255 * t288;
	t323 = t255 * t290;
	t240 = t301 ^ 2;
	t209 = t240 * t213 + 0.1e1;
	t315 = 0.2e1 * (-t240 * t339 - t301 * t329) / t209 ^ 2;
	t314 = -0.2e1 * t335;
	t312 = t223 * t334;
	t310 = -0.2e1 * t241 * t326;
	t309 = qJD(6) * t323 - t307;
	t304 = t330 * t221 - t332;
	t243 = t253 * t290 - t311;
	t303 = -t243 * t257 + t261 * t325;
	t298 = -qJD(4) * t324 + qJD(6) * t306 + t248 * t290;
	t264 = t284 * t271;
	t246 = t283 * t270 - t286 * t299;
	t238 = -qJD(4) * t260 - t265 * t290;
	t230 = t279 * t306 + t280 * t323;
	t229 = t279 * t323 - t280 * t306;
	t226 = -qJD(4) * t241 + t247 * t290;
	t216 = 0.1e1 / t218;
	t207 = 0.1e1 / t209;
	t205 = t233 * t337;
	t204 = t303 * t233;
	t201 = (-t231 * t252 - t232 * t267) * t288 + t308 * t205;
	t200 = t204 * t308 - t231 * t243 + t232 * t261;
	t198 = t303 * t336 + (t261 * t310 - t226 * t257 + (t225 * t261 + t237 * t243 + t238 * t241) * t258) * t233;
	t196 = t336 * t337 + (t302 * t317 + (-t267 * t310 - t246 * t257 + (-t225 * t267 + t237 * t252 - t241 * t264) * t258) * t288) * t233;
	t1 = [0, t196, 0, t198, 0, 0; 0, (-t201 * t333 - t212 * t324) * t315 + ((t248 * t288 + t255 * t317) * t212 + (-t329 + t338) * t201 + (-t324 * t199 + (-t267 * t317 - t196 * t241 - t205 * t225 - t264 * t288 + (-t205 * t260 - t252 * t288) * t203) * t327 + (-t252 * t317 - t196 * t260 - t205 * t237 - t246 * t288 + (t205 * t241 + t267 * t288) * t203) * t328) * t213) * t207, 0, (-t200 * t333 - t212 * t245) * t315 + (t200 * t338 + t228 * t212 + (-t245 * t199 - t200 * t227 + (-t198 * t241 - t204 * t225 + t238 + (-t204 * t260 - t243) * t203) * t327 + (-t198 * t260 - t204 * t237 - t226 + (t204 * t241 - t261) * t203) * t328) * t213) * t207, 0, 0; 0, 0.2e1 * (-t220 * t229 + t230 * t331) * t335 + (0.2e1 * t230 * t312 + t309 * t220 * t280 + t298 * t332 + (t223 * t279 * t309 - t230 * t210 - t229 * t211 - t298 * t330) * t221) * t216, 0, -t304 * t301 * t314 + (t304 * t227 - ((-qJD(6) * t220 - 0.2e1 * t312) * t280 + (t210 * t280 + (t211 - t316) * t279) * t221) * t301) * t216, 0, t314 + 0.2e1 * (t210 * t221 * t216 + (-t216 * t334 - t221 * t335) * t223) * t223;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end