% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0, 0; t9, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0, 0; -t8, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:30
	% EndTime: 2019-10-10 13:29:30
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
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:31
	% DurationCPUTime: 0.12s
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
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:31
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (136->42), mult. (408->90), div. (0->0), fcn. (571->12), ass. (0->50)
	t184 = sin(qJ(2));
	t185 = sin(qJ(1));
	t188 = cos(qJ(2));
	t189 = cos(qJ(1));
	t208 = cos(pkin(6));
	t192 = t189 * t208;
	t173 = t184 * t192 + t185 * t188;
	t183 = sin(qJ(3));
	t187 = cos(qJ(3));
	t172 = t185 * t184 - t188 * t192;
	t179 = sin(pkin(7));
	t181 = cos(pkin(7));
	t180 = sin(pkin(6));
	t202 = t180 * t189;
	t190 = t172 * t181 + t179 * t202;
	t159 = -t173 * t187 + t190 * t183;
	t166 = -t172 * t179 + t181 * t202;
	t182 = sin(qJ(4));
	t186 = cos(qJ(4));
	t212 = t159 * t182 - t166 * t186;
	t211 = t159 * t186 + t166 * t182;
	t206 = t179 * t180;
	t205 = t179 * t182;
	t204 = t179 * t186;
	t203 = t180 * t185;
	t201 = t181 * t183;
	t200 = t181 * t187;
	t199 = t183 * t184;
	t198 = t183 * t188;
	t197 = t184 * t187;
	t196 = t187 * t188;
	t195 = t184 * t206;
	t194 = t179 * t203;
	t193 = t185 * t208;
	t191 = t208 * t179;
	t157 = -t173 * t183 - t190 * t187;
	t175 = -t184 * t193 + t189 * t188;
	t174 = -t189 * t184 - t188 * t193;
	t171 = t208 * t181 - t188 * t206;
	t170 = (-t181 * t199 + t196) * t180;
	t168 = -t174 * t179 + t181 * t203;
	t165 = t183 * t191 + (t181 * t198 + t197) * t180;
	t164 = t187 * t191 + (t181 * t196 - t199) * t180;
	t163 = t174 * t187 - t175 * t201;
	t162 = -t172 * t187 - t173 * t201;
	t161 = t175 * t187 + (t174 * t181 + t194) * t183;
	t160 = -t174 * t200 + t175 * t183 - t187 * t194;
	t156 = t161 * t186 + t168 * t182;
	t155 = -t161 * t182 + t168 * t186;
	t1 = [t211, t163 * t186 + t175 * t205, -t160 * t186, t155, 0, 0; t156, t162 * t186 + t173 * t205, t157 * t186, t212, 0, 0; 0, t170 * t186 + t182 * t195, t164 * t186, -t165 * t182 + t171 * t186, 0, 0; -t212, -t163 * t182 + t175 * t204, t160 * t182, -t156, 0, 0; t155, -t162 * t182 + t173 * t204, -t157 * t182, t211, 0, 0; 0, -t170 * t182 + t186 * t195, -t164 * t182, -t165 * t186 - t171 * t182, 0, 0; t157, t174 * t183 + t175 * t200, t161, 0, 0, 0; t160, -t172 * t183 + t173 * t200, -t159, 0, 0, 0; 0, (t181 * t197 + t198) * t180, t165, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:31
	% EndTime: 2019-10-10 13:29:31
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (222->44), mult. (512->90), div. (0->0), fcn. (717->12), ass. (0->53)
	t208 = sin(qJ(2));
	t209 = sin(qJ(1));
	t211 = cos(qJ(2));
	t212 = cos(qJ(1));
	t231 = cos(pkin(6));
	t215 = t212 * t231;
	t195 = t208 * t215 + t209 * t211;
	t207 = sin(qJ(3));
	t210 = cos(qJ(3));
	t194 = t209 * t208 - t211 * t215;
	t204 = sin(pkin(7));
	t206 = cos(pkin(7));
	t205 = sin(pkin(6));
	t225 = t205 * t212;
	t213 = t194 * t206 + t204 * t225;
	t181 = -t195 * t210 + t213 * t207;
	t188 = -t194 * t204 + t206 * t225;
	t203 = qJ(4) + qJ(5);
	t201 = sin(t203);
	t202 = cos(t203);
	t173 = t181 * t201 - t188 * t202;
	t174 = t181 * t202 + t188 * t201;
	t229 = t201 * t204;
	t228 = t202 * t204;
	t227 = t204 * t205;
	t226 = t205 * t209;
	t224 = t206 * t207;
	t223 = t206 * t210;
	t222 = t207 * t208;
	t221 = t207 * t211;
	t220 = t208 * t210;
	t219 = t210 * t211;
	t218 = t208 * t227;
	t217 = t204 * t226;
	t216 = t209 * t231;
	t214 = t231 * t204;
	t179 = -t195 * t207 - t213 * t210;
	t197 = -t208 * t216 + t212 * t211;
	t196 = -t212 * t208 - t211 * t216;
	t193 = t231 * t206 - t211 * t227;
	t192 = (-t206 * t222 + t219) * t205;
	t190 = -t196 * t204 + t206 * t226;
	t187 = t207 * t214 + (t206 * t221 + t220) * t205;
	t186 = t210 * t214 + (t206 * t219 - t222) * t205;
	t185 = t196 * t210 - t197 * t224;
	t184 = -t194 * t210 - t195 * t224;
	t183 = t197 * t210 + (t196 * t206 + t217) * t207;
	t182 = -t196 * t223 + t197 * t207 - t210 * t217;
	t178 = -t187 * t202 - t193 * t201;
	t177 = -t187 * t201 + t193 * t202;
	t176 = t183 * t202 + t190 * t201;
	t175 = -t183 * t201 + t190 * t202;
	t1 = [t174, t185 * t202 + t197 * t229, -t182 * t202, t175, t175, 0; t176, t184 * t202 + t195 * t229, t179 * t202, t173, t173, 0; 0, t192 * t202 + t201 * t218, t186 * t202, t177, t177, 0; -t173, -t185 * t201 + t197 * t228, t182 * t201, -t176, -t176, 0; t175, -t184 * t201 + t195 * t228, -t179 * t201, t174, t174, 0; 0, -t192 * t201 + t202 * t218, -t186 * t201, t178, t178, 0; t179, t196 * t207 + t197 * t223, t183, 0, 0, 0; t182, -t194 * t207 + t195 * t223, -t181, 0, 0, 0; 0, (t206 * t220 + t221) * t205, t187, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:29:32
	% EndTime: 2019-10-10 13:29:33
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (437->71), mult. (1029->135), div. (0->0), fcn. (1421->14), ass. (0->76)
	t292 = sin(qJ(2));
	t295 = cos(qJ(2));
	t296 = cos(qJ(1));
	t323 = cos(pkin(6));
	t302 = t296 * t323;
	t324 = sin(qJ(1));
	t276 = t292 * t302 + t324 * t295;
	t291 = sin(qJ(3));
	t294 = cos(qJ(3));
	t275 = t324 * t292 - t295 * t302;
	t289 = cos(pkin(7));
	t287 = sin(pkin(7));
	t288 = sin(pkin(6));
	t312 = t288 * t296;
	t304 = t287 * t312;
	t298 = t275 * t289 + t304;
	t257 = -t276 * t294 + t298 * t291;
	t268 = -t275 * t287 + t289 * t312;
	t286 = qJ(4) + qJ(5);
	t284 = sin(t286);
	t285 = cos(t286);
	t246 = t257 * t285 + t268 * t284;
	t290 = sin(qJ(6));
	t328 = t246 * t290;
	t293 = cos(qJ(6));
	t327 = t246 * t293;
	t244 = t257 * t284 - t268 * t285;
	t322 = t244 * t290;
	t299 = t323 * t324;
	t277 = -t296 * t292 - t295 * t299;
	t278 = -t292 * t299 + t296 * t295;
	t303 = t288 * t324;
	t300 = t287 * t303;
	t259 = t278 * t294 + (t277 * t289 + t300) * t291;
	t297 = -t277 * t287 + t289 * t303;
	t247 = t259 * t284 - t297 * t285;
	t321 = t247 * t290;
	t301 = t323 * t287;
	t307 = t292 * t294;
	t308 = t291 * t295;
	t267 = t291 * t301 + (t289 * t308 + t307) * t288;
	t313 = t287 * t288;
	t274 = t323 * t289 - t295 * t313;
	t252 = -t267 * t284 + t274 * t285;
	t320 = t252 * t290;
	t319 = t276 * t291;
	t317 = t284 * t287;
	t316 = t285 * t287;
	t315 = t285 * t290;
	t314 = t285 * t293;
	t311 = t289 * t291;
	t310 = t289 * t294;
	t309 = t291 * t292;
	t306 = t294 * t295;
	t305 = t292 * t313;
	t273 = (-t289 * t309 + t306) * t288;
	t272 = (t289 * t307 + t308) * t288;
	t266 = -t294 * t301 + (-t289 * t306 + t309) * t288;
	t264 = t273 * t285 + t284 * t305;
	t263 = t277 * t294 - t278 * t311;
	t262 = t277 * t291 + t278 * t310;
	t261 = -t275 * t294 - t276 * t311;
	t260 = -t275 * t291 + t276 * t310;
	t258 = -t277 * t310 + t278 * t291 - t294 * t300;
	t256 = -t298 * t294 - t319;
	t254 = t275 * t310 + t294 * t304 + t319;
	t253 = t267 * t285 + t274 * t284;
	t251 = t252 * t293;
	t250 = t263 * t285 + t278 * t317;
	t249 = t261 * t285 + t276 * t317;
	t248 = t259 * t285 + t297 * t284;
	t243 = t247 * t293;
	t242 = t244 * t293;
	t241 = t248 * t293 + t258 * t290;
	t240 = -t248 * t290 + t258 * t293;
	t1 = [t256 * t290 + t327, t250 * t293 + t262 * t290, -t258 * t314 + t259 * t290, -t243, -t243, t240; t241, t249 * t293 + t260 * t290, -t254 * t314 - t257 * t290, t242, t242, t254 * t293 + t328; 0, t264 * t293 + t272 * t290, -t266 * t314 + t267 * t290, t251, t251, -t253 * t290 + t266 * t293; t256 * t293 - t328, -t250 * t290 + t262 * t293, t258 * t315 + t259 * t293, t321, t321, -t241; t240, -t249 * t290 + t260 * t293, t254 * t315 - t257 * t293, -t322, -t322, -t254 * t290 + t327; 0, -t264 * t290 + t272 * t293, t266 * t315 + t267 * t293, -t320, -t320, -t253 * t293 - t266 * t290; t244, t263 * t284 - t278 * t316, -t258 * t284, t248, t248, 0; t247, t261 * t284 - t276 * t316, -t254 * t284, -t246, -t246, 0; 0, t273 * t284 - t285 * t305, -t266 * t284, t253, t253, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end