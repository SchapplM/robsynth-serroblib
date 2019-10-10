% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:22
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_jacobiR_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:20
	% EndTime: 2019-10-09 21:22:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (20->14), mult. (62->33), div. (0->0), fcn. (88->10), ass. (0->20)
	t53 = sin(pkin(13));
	t59 = cos(pkin(6));
	t68 = t53 * t59;
	t54 = sin(pkin(7));
	t55 = sin(pkin(6));
	t67 = t54 * t55;
	t66 = t54 * t59;
	t56 = cos(pkin(14));
	t58 = cos(pkin(7));
	t65 = t56 * t58;
	t57 = cos(pkin(13));
	t64 = t57 * t59;
	t52 = sin(pkin(14));
	t63 = -(-t53 * t52 + t56 * t64) * t58 + t57 * t67;
	t62 = (-t57 * t52 - t56 * t68) * t58 + t53 * t67;
	t61 = cos(qJ(3));
	t60 = sin(qJ(3));
	t51 = -t52 * t68 + t57 * t56;
	t49 = t52 * t64 + t53 * t56;
	t1 = [0, 0, -t51 * t60 + t62 * t61, 0, 0, 0; 0, 0, -t49 * t60 - t63 * t61, 0, 0, 0; 0, 0, t61 * t66 + (-t52 * t60 + t61 * t65) * t55, 0, 0, 0; 0, 0, -t51 * t61 - t62 * t60, 0, 0, 0; 0, 0, -t49 * t61 + t63 * t60, 0, 0, 0; 0, 0, -t60 * t66 + (-t52 * t61 - t60 * t65) * t55, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:21
	% DurationCPUTime: 0.18s
	% Computational Cost: add. (121->33), mult. (360->77), div. (0->0), fcn. (493->14), ass. (0->38)
	t132 = sin(pkin(13));
	t140 = cos(pkin(6));
	t157 = t132 * t140;
	t134 = sin(pkin(7));
	t135 = sin(pkin(6));
	t156 = t134 * t135;
	t155 = t134 * t140;
	t139 = cos(pkin(7));
	t154 = t135 * t139;
	t137 = cos(pkin(13));
	t153 = t137 * t140;
	t138 = cos(pkin(8));
	t141 = sin(qJ(4));
	t152 = t138 * t141;
	t143 = cos(qJ(4));
	t151 = t138 * t143;
	t142 = sin(qJ(3));
	t150 = t139 * t142;
	t149 = t137 * t156;
	t131 = sin(pkin(14));
	t136 = cos(pkin(14));
	t126 = -t131 * t132 + t136 * t153;
	t127 = t131 * t153 + t132 * t136;
	t144 = cos(qJ(3));
	t117 = -t127 * t142 + (t126 * t139 - t149) * t144;
	t133 = sin(pkin(8));
	t148 = t117 * t138 + (-t126 * t134 - t137 * t154) * t133;
	t129 = -t131 * t157 + t136 * t137;
	t128 = -t131 * t137 - t136 * t157;
	t145 = t128 * t139 + t132 * t156;
	t119 = -t129 * t142 + t144 * t145;
	t147 = t119 * t138 + (-t128 * t134 + t132 * t154) * t133;
	t121 = t144 * t155 + (t136 * t139 * t144 - t131 * t142) * t135;
	t146 = t121 * t138 + (-t136 * t156 + t139 * t140) * t133;
	t122 = t142 * t155 + (t131 * t144 + t136 * t150) * t135;
	t120 = t129 * t144 + t142 * t145;
	t118 = t126 * t150 + t127 * t144 - t142 * t149;
	t1 = [0, 0, t119 * t143 - t120 * t152, -t120 * t141 + t143 * t147, 0, 0; 0, 0, t117 * t143 - t118 * t152, -t118 * t141 + t143 * t148, 0, 0; 0, 0, t121 * t143 - t122 * t152, -t122 * t141 + t143 * t146, 0, 0; 0, 0, -t119 * t141 - t120 * t151, -t120 * t143 - t141 * t147, 0, 0; 0, 0, -t117 * t141 - t118 * t151, -t118 * t143 - t141 * t148, 0, 0; 0, 0, -t121 * t141 - t122 * t151, -t122 * t143 - t141 * t146, 0, 0; 0, 0, t120 * t133, 0, 0, 0; 0, 0, t118 * t133, 0, 0, 0; 0, 0, t122 * t133, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:21
	% EndTime: 2019-10-09 21:22:22
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (295->51), mult. (875->112), div. (0->0), fcn. (1191->16), ass. (0->57)
	t205 = sin(pkin(13));
	t213 = cos(pkin(6));
	t234 = t205 * t213;
	t206 = sin(pkin(8));
	t214 = sin(qJ(5));
	t233 = t206 * t214;
	t217 = cos(qJ(5));
	t232 = t206 * t217;
	t207 = sin(pkin(7));
	t208 = sin(pkin(6));
	t231 = t207 * t208;
	t230 = t207 * t213;
	t212 = cos(pkin(7));
	t229 = t208 * t212;
	t210 = cos(pkin(13));
	t228 = t210 * t213;
	t211 = cos(pkin(8));
	t215 = sin(qJ(4));
	t227 = t211 * t215;
	t218 = cos(qJ(4));
	t226 = t211 * t218;
	t216 = sin(qJ(3));
	t225 = t212 * t216;
	t224 = t210 * t231;
	t204 = sin(pkin(14));
	t209 = cos(pkin(14));
	t199 = -t205 * t204 + t209 * t228;
	t200 = t204 * t228 + t205 * t209;
	t219 = cos(qJ(3));
	t189 = -t200 * t216 + (t199 * t212 - t224) * t219;
	t196 = -t199 * t207 - t210 * t229;
	t223 = t189 * t211 + t196 * t206;
	t202 = -t204 * t234 + t210 * t209;
	t201 = -t210 * t204 - t209 * t234;
	t220 = t201 * t212 + t205 * t231;
	t191 = -t202 * t216 + t220 * t219;
	t197 = -t201 * t207 + t205 * t229;
	t222 = t191 * t211 + t197 * t206;
	t194 = t219 * t230 + (t209 * t212 * t219 - t204 * t216) * t208;
	t198 = -t209 * t231 + t213 * t212;
	t221 = t194 * t211 + t198 * t206;
	t195 = t216 * t230 + (t204 * t219 + t209 * t225) * t208;
	t193 = -t194 * t206 + t198 * t211;
	t192 = t202 * t219 + t220 * t216;
	t190 = t199 * t225 + t200 * t219 - t216 * t224;
	t188 = t194 * t218 - t195 * t227;
	t187 = -t191 * t206 + t197 * t211;
	t186 = -t189 * t206 + t196 * t211;
	t185 = t195 * t218 + t221 * t215;
	t184 = -t195 * t215 + t221 * t218;
	t183 = t191 * t218 - t192 * t227;
	t182 = t189 * t218 - t190 * t227;
	t181 = t192 * t218 + t222 * t215;
	t180 = -t192 * t215 + t222 * t218;
	t179 = t190 * t218 + t223 * t215;
	t178 = -t190 * t215 + t223 * t218;
	t1 = [0, 0, t183 * t217 + t192 * t233, t180 * t217, -t181 * t214 + t187 * t217, 0; 0, 0, t182 * t217 + t190 * t233, t178 * t217, -t179 * t214 + t186 * t217, 0; 0, 0, t188 * t217 + t195 * t233, t184 * t217, -t185 * t214 + t193 * t217, 0; 0, 0, -t183 * t214 + t192 * t232, -t180 * t214, -t181 * t217 - t187 * t214, 0; 0, 0, -t182 * t214 + t190 * t232, -t178 * t214, -t179 * t217 - t186 * t214, 0; 0, 0, -t188 * t214 + t195 * t232, -t184 * t214, -t185 * t217 - t193 * t214, 0; 0, 0, t191 * t215 + t192 * t226, t181, 0, 0; 0, 0, t189 * t215 + t190 * t226, t179, 0, 0; 0, 0, t194 * t215 + t195 * t226, t185, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:22:22
	% EndTime: 2019-10-09 21:22:22
	% DurationCPUTime: 0.37s
	% Computational Cost: add. (652->75), mult. (1920->156), div. (0->0), fcn. (2602->18), ass. (0->73)
	t274 = sin(pkin(13));
	t282 = cos(pkin(6));
	t310 = t274 * t282;
	t275 = sin(pkin(8));
	t284 = sin(qJ(5));
	t309 = t275 * t284;
	t288 = cos(qJ(5));
	t308 = t275 * t288;
	t276 = sin(pkin(7));
	t277 = sin(pkin(6));
	t307 = t276 * t277;
	t306 = t276 * t282;
	t281 = cos(pkin(7));
	t305 = t277 * t281;
	t279 = cos(pkin(13));
	t304 = t279 * t282;
	t280 = cos(pkin(8));
	t285 = sin(qJ(4));
	t303 = t280 * t285;
	t289 = cos(qJ(4));
	t302 = t280 * t289;
	t286 = sin(qJ(3));
	t301 = t281 * t286;
	t283 = sin(qJ(6));
	t300 = t283 * t288;
	t287 = cos(qJ(6));
	t299 = t287 * t288;
	t298 = t279 * t307;
	t273 = sin(pkin(14));
	t278 = cos(pkin(14));
	t268 = -t273 * t274 + t278 * t304;
	t297 = -t268 * t276 - t279 * t305;
	t270 = -t273 * t279 - t278 * t310;
	t296 = -t270 * t276 + t274 * t305;
	t295 = t270 * t281 + t274 * t307;
	t294 = -t278 * t307 + t281 * t282;
	t293 = t297 * t275;
	t292 = t296 * t275;
	t291 = t294 * t275;
	t290 = cos(qJ(3));
	t271 = -t273 * t310 + t278 * t279;
	t269 = t273 * t304 + t274 * t278;
	t266 = t286 * t306 + (t273 * t290 + t278 * t301) * t277;
	t265 = t290 * t306 + (t278 * t281 * t290 - t273 * t286) * t277;
	t261 = -t265 * t275 + t294 * t280;
	t260 = t271 * t290 + t295 * t286;
	t259 = -t271 * t286 + t295 * t290;
	t258 = t268 * t301 + t269 * t290 - t286 * t298;
	t257 = -t269 * t286 + (t268 * t281 - t298) * t290;
	t254 = t265 * t289 - t266 * t303;
	t253 = t265 * t285 + t266 * t302;
	t252 = -t259 * t275 + t296 * t280;
	t251 = -t257 * t275 + t297 * t280;
	t250 = t266 * t289 + (t265 * t280 + t291) * t285;
	t249 = -t265 * t302 + t266 * t285 - t289 * t291;
	t248 = t254 * t288 + t266 * t309;
	t247 = t259 * t289 - t260 * t303;
	t246 = t259 * t285 + t260 * t302;
	t245 = t257 * t289 - t258 * t303;
	t244 = t257 * t285 + t258 * t302;
	t243 = t250 * t288 + t261 * t284;
	t242 = -t250 * t284 + t261 * t288;
	t241 = t260 * t289 + (t259 * t280 + t292) * t285;
	t240 = -t259 * t302 + t260 * t285 - t289 * t292;
	t239 = t258 * t289 + (t257 * t280 + t293) * t285;
	t238 = -t257 * t302 + t258 * t285 - t289 * t293;
	t237 = t247 * t288 + t260 * t309;
	t236 = t245 * t288 + t258 * t309;
	t235 = t241 * t288 + t252 * t284;
	t234 = -t241 * t284 + t252 * t288;
	t233 = t239 * t288 + t251 * t284;
	t232 = -t239 * t284 + t251 * t288;
	t1 = [0, 0, t237 * t287 + t246 * t283, -t240 * t299 + t241 * t283, t234 * t287, -t235 * t283 + t240 * t287; 0, 0, t236 * t287 + t244 * t283, -t238 * t299 + t239 * t283, t232 * t287, -t233 * t283 + t238 * t287; 0, 0, t248 * t287 + t253 * t283, -t249 * t299 + t250 * t283, t242 * t287, -t243 * t283 + t249 * t287; 0, 0, -t237 * t283 + t246 * t287, t240 * t300 + t241 * t287, -t234 * t283, -t235 * t287 - t240 * t283; 0, 0, -t236 * t283 + t244 * t287, t238 * t300 + t239 * t287, -t232 * t283, -t233 * t287 - t238 * t283; 0, 0, -t248 * t283 + t253 * t287, t249 * t300 + t250 * t287, -t242 * t283, -t243 * t287 - t249 * t283; 0, 0, t247 * t284 - t260 * t308, -t240 * t284, t235, 0; 0, 0, t245 * t284 - t258 * t308, -t238 * t284, t233, 0; 0, 0, t254 * t284 - t266 * t308, -t249 * t284, t243, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end