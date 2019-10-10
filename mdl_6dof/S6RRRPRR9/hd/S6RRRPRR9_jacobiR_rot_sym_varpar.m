% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:08
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR9_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR9_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t50 = sin(qJ(2));
	t51 = sin(qJ(1));
	t57 = t51 * t50;
	t52 = cos(qJ(2));
	t56 = t51 * t52;
	t53 = cos(qJ(1));
	t55 = t53 * t50;
	t54 = t53 * t52;
	t49 = cos(pkin(6));
	t48 = sin(pkin(6));
	t47 = -t49 * t57 + t54;
	t46 = -t49 * t56 - t55;
	t45 = -t49 * t55 - t56;
	t44 = -t49 * t54 + t57;
	t1 = [t45, t46, 0, 0, 0, 0; t47, -t44, 0, 0, 0, 0; 0, t48 * t52, 0, 0, 0, 0; t44, -t47, 0, 0, 0, 0; t46, t45, 0, 0, 0, 0; 0, -t48 * t50, 0, 0, 0, 0; t53 * t48, 0, 0, 0, 0, 0; t51 * t48, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:08
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (58->25), mult. (178->58), div. (0->0), fcn. (255->10), ass. (0->32)
	t115 = sin(qJ(2));
	t116 = sin(qJ(1));
	t118 = cos(qJ(2));
	t119 = cos(qJ(1));
	t136 = cos(pkin(6));
	t124 = t119 * t136;
	t104 = t116 * t115 - t118 * t124;
	t105 = t115 * t124 + t116 * t118;
	t113 = cos(pkin(7));
	t114 = sin(qJ(3));
	t117 = cos(qJ(3));
	t111 = sin(pkin(7));
	t112 = sin(pkin(6));
	t133 = t112 * t119;
	t126 = t111 * t133;
	t137 = (t104 * t113 + t126) * t117 + t105 * t114;
	t134 = t112 * t116;
	t132 = t113 * t114;
	t131 = t113 * t117;
	t130 = t114 * t115;
	t129 = t114 * t118;
	t128 = t115 * t117;
	t127 = t117 * t118;
	t125 = t116 * t136;
	t123 = t136 * t111;
	t106 = -t119 * t115 - t118 * t125;
	t121 = t106 * t113 + t111 * t134;
	t120 = t104 * t132 - t105 * t117 + t114 * t126;
	t107 = -t115 * t125 + t119 * t118;
	t103 = t107 * t117 + t121 * t114;
	t102 = -t107 * t114 + t121 * t117;
	t1 = [t120, t106 * t117 - t107 * t132, t102, 0, 0, 0; t103, -t104 * t117 - t105 * t132, -t137, 0, 0, 0; 0, (-t113 * t130 + t127) * t112, t117 * t123 + (t113 * t127 - t130) * t112, 0, 0, 0; t137, -t106 * t114 - t107 * t131, -t103, 0, 0, 0; t102, t104 * t114 - t105 * t131, t120, 0, 0, 0; 0, (-t113 * t128 - t129) * t112, -t114 * t123 + (-t113 * t129 - t128) * t112, 0, 0, 0; -t104 * t111 + t113 * t133, t107 * t111, 0, 0, 0, 0; -t106 * t111 + t113 * t134, t105 * t111, 0, 0, 0, 0; 0, t112 * t115 * t111, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:08
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (100->28), mult. (288->62), div. (0->0), fcn. (407->12), ass. (0->31)
	t146 = cos(pkin(6));
	t129 = sin(pkin(6));
	t134 = sin(qJ(1));
	t145 = t129 * t134;
	t137 = cos(qJ(1));
	t144 = t129 * t137;
	t143 = t134 * t146;
	t142 = t137 * t146;
	t127 = sin(pkin(13));
	t130 = cos(pkin(13));
	t132 = sin(qJ(3));
	t135 = cos(qJ(3));
	t141 = t135 * t127 + t132 * t130;
	t122 = t132 * t127 - t135 * t130;
	t128 = sin(pkin(7));
	t114 = t122 * t128;
	t131 = cos(pkin(7));
	t116 = t122 * t131;
	t133 = sin(qJ(2));
	t136 = cos(qJ(2));
	t118 = t134 * t133 - t136 * t142;
	t119 = t133 * t142 + t134 * t136;
	t140 = -t114 * t144 - t118 * t116 + t119 * t141;
	t115 = t141 * t128;
	t117 = t141 * t131;
	t139 = t115 * t144 + t118 * t117 + t119 * t122;
	t120 = -t137 * t133 - t136 * t143;
	t121 = -t133 * t143 + t137 * t136;
	t138 = t115 * t145 + t120 * t117 - t121 * t122;
	t113 = -t114 * t145 - t120 * t116 - t121 * t141;
	t1 = [t139, -t121 * t117 - t120 * t122, t113, 0, 0, 0; t138, -t119 * t117 + t118 * t122, -t140, 0, 0, 0; 0, (-t117 * t133 - t122 * t136) * t129, -t146 * t114 + (-t116 * t136 - t133 * t141) * t129, 0, 0, 0; t140, t121 * t116 - t120 * t141, -t138, 0, 0, 0; t113, t119 * t116 + t118 * t141, t139, 0, 0, 0; 0, (t116 * t133 - t136 * t141) * t129, -t146 * t115 + (-t117 * t136 + t122 * t133) * t129, 0, 0, 0; -t118 * t128 + t131 * t144, t121 * t128, 0, 0, 0, 0; -t120 * t128 + t131 * t145, t119 * t128, 0, 0, 0, 0; 0, t129 * t133 * t128, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:08
	% EndTime: 2019-10-10 12:08:08
	% DurationCPUTime: 0.32s
	% Computational Cost: add. (217->46), mult. (621->94), div. (0->0), fcn. (865->14), ass. (0->49)
	t198 = sin(pkin(7));
	t197 = sin(pkin(13));
	t202 = sin(qJ(3));
	t218 = cos(pkin(13));
	t220 = cos(qJ(3));
	t209 = t220 * t197 + t202 * t218;
	t184 = t209 * t198;
	t200 = cos(pkin(7));
	t186 = t209 * t200;
	t203 = sin(qJ(2));
	t204 = sin(qJ(1));
	t206 = cos(qJ(2));
	t207 = cos(qJ(1));
	t219 = cos(pkin(6));
	t210 = t207 * t219;
	t188 = t204 * t203 - t206 * t210;
	t189 = t203 * t210 + t204 * t206;
	t208 = -t202 * t197 + t220 * t218;
	t199 = sin(pkin(6));
	t213 = t199 * t207;
	t172 = t184 * t213 + t188 * t186 - t189 * t208;
	t180 = -t188 * t198 + t200 * t213;
	t201 = sin(qJ(5));
	t205 = cos(qJ(5));
	t222 = t172 * t205 + t180 * t201;
	t221 = t172 * t201 - t180 * t205;
	t217 = t198 * t199;
	t216 = t198 * t201;
	t215 = t198 * t205;
	t214 = t199 * t204;
	t212 = t203 * t217;
	t211 = t204 * t219;
	t183 = t208 * t198;
	t185 = t208 * t200;
	t170 = -t183 * t213 - t188 * t185 - t189 * t209;
	t190 = -t207 * t203 - t206 * t211;
	t191 = -t203 * t211 + t207 * t206;
	t174 = t184 * t214 + t190 * t186 + t191 * t208;
	t176 = t219 * t184 + (t186 * t206 + t203 * t208) * t199;
	t187 = t219 * t200 - t206 * t217;
	t182 = -t190 * t198 + t200 * t214;
	t179 = (-t186 * t203 + t206 * t208) * t199;
	t178 = -t191 * t186 + t190 * t208;
	t177 = -t189 * t186 - t188 * t208;
	t175 = t219 * t183 + (t185 * t206 - t203 * t209) * t199;
	t173 = t183 * t214 + t190 * t185 - t191 * t209;
	t169 = t174 * t205 + t182 * t201;
	t168 = -t174 * t201 + t182 * t205;
	t1 = [t222, t178 * t205 + t191 * t216, t173 * t205, 0, t168, 0; t169, t177 * t205 + t189 * t216, t170 * t205, 0, t221, 0; 0, t179 * t205 + t201 * t212, t175 * t205, 0, -t176 * t201 + t187 * t205, 0; -t221, -t178 * t201 + t191 * t215, -t173 * t201, 0, -t169, 0; t168, -t177 * t201 + t189 * t215, -t170 * t201, 0, t222, 0; 0, -t179 * t201 + t205 * t212, -t175 * t201, 0, -t176 * t205 - t187 * t201, 0; t170, t191 * t185 + t190 * t209, t174, 0, 0, 0; -t173, t189 * t185 - t188 * t209, -t172, 0, 0, 0; 0, (t185 * t203 + t206 * t209) * t199, t176, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:09
	% EndTime: 2019-10-10 12:08:09
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (461->63), mult. (1310->135), div. (0->0), fcn. (1807->16), ass. (0->65)
	t270 = sin(pkin(13));
	t276 = sin(qJ(3));
	t298 = cos(pkin(13));
	t300 = cos(qJ(3));
	t266 = -t300 * t270 - t276 * t298;
	t271 = sin(pkin(7));
	t257 = t266 * t271;
	t273 = cos(pkin(7));
	t259 = t266 * t273;
	t277 = sin(qJ(2));
	t278 = sin(qJ(1));
	t281 = cos(qJ(2));
	t282 = cos(qJ(1));
	t299 = cos(pkin(6));
	t288 = t282 * t299;
	t261 = t278 * t277 - t281 * t288;
	t262 = t277 * t288 + t278 * t281;
	t285 = -t276 * t270 + t300 * t298;
	t272 = sin(pkin(6));
	t293 = t272 * t282;
	t238 = -t257 * t293 - t261 * t259 - t262 * t285;
	t253 = -t261 * t271 + t273 * t293;
	t275 = sin(qJ(5));
	t280 = cos(qJ(5));
	t227 = t238 * t280 + t253 * t275;
	t274 = sin(qJ(6));
	t279 = cos(qJ(6));
	t256 = t285 * t271;
	t258 = t285 * t273;
	t286 = -t256 * t293 - t261 * t258 + t262 * t266;
	t304 = t227 * t274 - t279 * t286;
	t303 = t227 * t279 + t274 * t286;
	t225 = t238 * t275 - t253 * t280;
	t297 = t271 * t272;
	t296 = t271 * t275;
	t295 = t271 * t280;
	t294 = t272 * t278;
	t292 = t274 * t280;
	t291 = t279 * t280;
	t290 = t277 * t297;
	t289 = t278 * t299;
	t263 = -t282 * t277 - t281 * t289;
	t287 = -t263 * t271 + t273 * t294;
	t264 = -t277 * t289 + t282 * t281;
	t284 = -t257 * t294 - t263 * t259 + t264 * t285;
	t283 = -t299 * t257 + (-t259 * t281 + t277 * t285) * t272;
	t260 = t299 * t273 - t281 * t297;
	t251 = (t259 * t277 + t281 * t285) * t272;
	t250 = (t258 * t277 - t266 * t281) * t272;
	t249 = t251 * t280 + t275 * t290;
	t248 = t264 * t259 + t263 * t285;
	t247 = t264 * t258 - t263 * t266;
	t246 = t262 * t259 - t261 * t285;
	t245 = t262 * t258 + t261 * t266;
	t243 = t299 * t256 + (t258 * t281 + t266 * t277) * t272;
	t240 = t256 * t294 + t263 * t258 + t264 * t266;
	t233 = t248 * t280 + t264 * t296;
	t232 = t246 * t280 + t262 * t296;
	t231 = t260 * t275 + t280 * t283;
	t230 = t260 * t280 - t275 * t283;
	t229 = t287 * t275 + t280 * t284;
	t228 = t275 * t284 - t287 * t280;
	t224 = t229 * t279 - t240 * t274;
	t223 = -t229 * t274 - t240 * t279;
	t1 = [t303, t233 * t279 + t247 * t274, t240 * t291 + t274 * t284, 0, -t228 * t279, t223; t224, t232 * t279 + t245 * t274, -t238 * t274 + t286 * t291, 0, t225 * t279, t304; 0, t249 * t279 + t250 * t274, t243 * t291 + t274 * t283, 0, t230 * t279, -t231 * t274 - t243 * t279; -t304, -t233 * t274 + t247 * t279, -t240 * t292 + t279 * t284, 0, t228 * t274, -t224; t223, -t232 * t274 + t245 * t279, -t238 * t279 - t286 * t292, 0, -t225 * t274, t303; 0, -t249 * t274 + t250 * t279, -t243 * t292 + t279 * t283, 0, -t230 * t274, -t231 * t279 + t243 * t274; t225, t248 * t275 - t264 * t295, t240 * t275, 0, t229, 0; t228, t246 * t275 - t262 * t295, t286 * t275, 0, -t227, 0; 0, t251 * t275 - t280 * t290, t243 * t275, 0, t231, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end