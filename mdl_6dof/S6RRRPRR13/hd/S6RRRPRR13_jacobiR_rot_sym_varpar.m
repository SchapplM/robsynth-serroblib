% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR13
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
% Datum: 2019-10-10 12:14
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR13_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR13_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
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
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
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
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:18
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
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:18
	% DurationCPUTime: 0.20s
	% Computational Cost: add. (100->37), mult. (304->83), div. (0->0), fcn. (425->12), ass. (0->43)
	t168 = sin(qJ(2));
	t169 = sin(qJ(1));
	t171 = cos(qJ(2));
	t172 = cos(qJ(1));
	t190 = cos(pkin(6));
	t174 = t172 * t190;
	t156 = t168 * t174 + t169 * t171;
	t167 = sin(qJ(3));
	t170 = cos(qJ(3));
	t155 = t168 * t169 - t171 * t174;
	t163 = sin(pkin(7));
	t166 = cos(pkin(7));
	t164 = sin(pkin(6));
	t185 = t164 * t172;
	t173 = t155 * t166 + t163 * t185;
	t145 = -t156 * t170 + t173 * t167;
	t162 = sin(pkin(13));
	t188 = t162 * t163;
	t165 = cos(pkin(13));
	t187 = t163 * t165;
	t186 = t164 * t169;
	t184 = t166 * t167;
	t183 = t166 * t170;
	t182 = t167 * t168;
	t181 = t167 * t171;
	t180 = t168 * t170;
	t179 = t170 * t171;
	t178 = t163 * t164 * t168;
	t177 = t163 * t186;
	t176 = t163 * t190;
	t175 = t169 * t190;
	t144 = -t156 * t167 - t173 * t170;
	t158 = -t168 * t175 + t171 * t172;
	t157 = -t172 * t168 - t171 * t175;
	t154 = (-t166 * t182 + t179) * t164;
	t152 = -t157 * t163 + t166 * t186;
	t151 = -t155 * t163 + t166 * t185;
	t150 = t170 * t176 + (t166 * t179 - t182) * t164;
	t149 = t157 * t170 - t158 * t184;
	t148 = -t155 * t170 - t156 * t184;
	t147 = t158 * t170 + (t157 * t166 + t177) * t167;
	t146 = -t157 * t183 + t158 * t167 - t170 * t177;
	t1 = [t145 * t165 + t151 * t162, t149 * t165 + t158 * t188, -t146 * t165, 0, 0, 0; t147 * t165 + t152 * t162, t148 * t165 + t156 * t188, t144 * t165, 0, 0, 0; 0, t154 * t165 + t162 * t178, t150 * t165, 0, 0, 0; -t145 * t162 + t151 * t165, -t149 * t162 + t158 * t187, t146 * t162, 0, 0, 0; -t147 * t162 + t152 * t165, -t148 * t162 + t156 * t187, -t144 * t162, 0, 0, 0; 0, -t154 * t162 + t165 * t178, -t150 * t162, 0, 0, 0; t144, t157 * t167 + t158 * t183, t147, 0, 0, 0; t146, -t155 * t167 + t156 * t183, -t145, 0, 0, 0; 0, (t166 * t180 + t181) * t164, t167 * t176 + (t166 * t181 + t180) * t164, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:18
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (174->43), mult. (408->90), div. (0->0), fcn. (571->12), ass. (0->51)
	t193 = sin(qJ(2));
	t194 = sin(qJ(1));
	t196 = cos(qJ(2));
	t197 = cos(qJ(1));
	t216 = cos(pkin(6));
	t200 = t197 * t216;
	t180 = t193 * t200 + t194 * t196;
	t192 = sin(qJ(3));
	t195 = cos(qJ(3));
	t179 = t194 * t193 - t196 * t200;
	t189 = sin(pkin(7));
	t191 = cos(pkin(7));
	t190 = sin(pkin(6));
	t210 = t190 * t197;
	t198 = t179 * t191 + t189 * t210;
	t166 = -t180 * t195 + t198 * t192;
	t173 = -t179 * t189 + t191 * t210;
	t188 = pkin(13) + qJ(5);
	t186 = sin(t188);
	t187 = cos(t188);
	t220 = t166 * t186 - t173 * t187;
	t219 = t166 * t187 + t173 * t186;
	t214 = t186 * t189;
	t213 = t187 * t189;
	t212 = t189 * t190;
	t211 = t190 * t194;
	t209 = t191 * t192;
	t208 = t191 * t195;
	t207 = t192 * t193;
	t206 = t192 * t196;
	t205 = t193 * t195;
	t204 = t195 * t196;
	t203 = t193 * t212;
	t202 = t189 * t211;
	t201 = t194 * t216;
	t199 = t216 * t189;
	t164 = -t180 * t192 - t198 * t195;
	t182 = -t193 * t201 + t197 * t196;
	t181 = -t197 * t193 - t196 * t201;
	t178 = t216 * t191 - t196 * t212;
	t177 = (-t191 * t207 + t204) * t190;
	t175 = -t181 * t189 + t191 * t211;
	t172 = t192 * t199 + (t191 * t206 + t205) * t190;
	t171 = t195 * t199 + (t191 * t204 - t207) * t190;
	t170 = t181 * t195 - t182 * t209;
	t169 = -t179 * t195 - t180 * t209;
	t168 = t182 * t195 + (t181 * t191 + t202) * t192;
	t167 = -t181 * t208 + t182 * t192 - t195 * t202;
	t163 = t168 * t187 + t175 * t186;
	t162 = -t168 * t186 + t175 * t187;
	t1 = [t219, t170 * t187 + t182 * t214, -t167 * t187, 0, t162, 0; t163, t169 * t187 + t180 * t214, t164 * t187, 0, t220, 0; 0, t177 * t187 + t186 * t203, t171 * t187, 0, -t172 * t186 + t178 * t187, 0; -t220, -t170 * t186 + t182 * t213, t167 * t186, 0, -t163, 0; t162, -t169 * t186 + t180 * t213, -t164 * t186, 0, t219, 0; 0, -t177 * t186 + t187 * t203, -t171 * t186, 0, -t172 * t187 - t178 * t186, 0; t164, t181 * t192 + t182 * t208, t168, 0, 0, 0; t167, -t179 * t192 + t180 * t208, -t166, 0, 0, 0; 0, (t191 * t205 + t206) * t190, t172, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:19
	% EndTime: 2019-10-10 12:14:19
	% DurationCPUTime: 0.58s
	% Computational Cost: add. (362->67), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->70)
	t267 = sin(qJ(2));
	t268 = sin(qJ(1));
	t271 = cos(qJ(2));
	t272 = cos(qJ(1));
	t296 = cos(pkin(6));
	t276 = t272 * t296;
	t251 = t267 * t276 + t268 * t271;
	t266 = sin(qJ(3));
	t270 = cos(qJ(3));
	t250 = t268 * t267 - t271 * t276;
	t264 = cos(pkin(7));
	t262 = sin(pkin(7));
	t263 = sin(pkin(6));
	t287 = t263 * t272;
	t278 = t262 * t287;
	t274 = t250 * t264 + t278;
	t232 = -t251 * t270 + t274 * t266;
	t243 = -t250 * t262 + t264 * t287;
	t261 = pkin(13) + qJ(5);
	t259 = sin(t261);
	t260 = cos(t261);
	t222 = t232 * t260 + t243 * t259;
	t265 = sin(qJ(6));
	t300 = t222 * t265;
	t269 = cos(qJ(6));
	t299 = t222 * t269;
	t220 = t232 * t259 - t243 * t260;
	t295 = t251 * t266;
	t293 = t259 * t262;
	t292 = t260 * t262;
	t291 = t260 * t265;
	t290 = t260 * t269;
	t289 = t262 * t263;
	t288 = t263 * t268;
	t286 = t264 * t266;
	t285 = t264 * t270;
	t284 = t266 * t267;
	t283 = t266 * t271;
	t282 = t267 * t270;
	t281 = t270 * t271;
	t280 = t267 * t289;
	t279 = t262 * t288;
	t277 = t268 * t296;
	t275 = t296 * t262;
	t252 = -t272 * t267 - t271 * t277;
	t273 = -t252 * t262 + t264 * t288;
	t253 = -t267 * t277 + t272 * t271;
	t249 = t296 * t264 - t271 * t289;
	t248 = (-t264 * t284 + t281) * t263;
	t247 = (t264 * t282 + t283) * t263;
	t242 = t266 * t275 + (t264 * t283 + t282) * t263;
	t241 = -t270 * t275 + (-t264 * t281 + t284) * t263;
	t239 = t252 * t270 - t253 * t286;
	t238 = t252 * t266 + t253 * t285;
	t237 = -t250 * t270 - t251 * t286;
	t236 = -t250 * t266 + t251 * t285;
	t235 = t248 * t260 + t259 * t280;
	t234 = t253 * t270 + (t252 * t264 + t279) * t266;
	t233 = -t252 * t285 + t253 * t266 - t270 * t279;
	t231 = -t274 * t270 - t295;
	t229 = t250 * t285 + t270 * t278 + t295;
	t228 = t242 * t260 + t249 * t259;
	t227 = -t242 * t259 + t249 * t260;
	t226 = t239 * t260 + t253 * t293;
	t225 = t237 * t260 + t251 * t293;
	t224 = t234 * t260 + t273 * t259;
	t223 = t234 * t259 - t273 * t260;
	t219 = t224 * t269 + t233 * t265;
	t218 = -t224 * t265 + t233 * t269;
	t1 = [t231 * t265 + t299, t226 * t269 + t238 * t265, -t233 * t290 + t234 * t265, 0, -t223 * t269, t218; t219, t225 * t269 + t236 * t265, -t229 * t290 - t232 * t265, 0, t220 * t269, t229 * t269 + t300; 0, t235 * t269 + t247 * t265, -t241 * t290 + t242 * t265, 0, t227 * t269, -t228 * t265 + t241 * t269; t231 * t269 - t300, -t226 * t265 + t238 * t269, t233 * t291 + t234 * t269, 0, t223 * t265, -t219; t218, -t225 * t265 + t236 * t269, t229 * t291 - t232 * t269, 0, -t220 * t265, -t229 * t265 + t299; 0, -t235 * t265 + t247 * t269, t241 * t291 + t242 * t269, 0, -t227 * t265, -t228 * t269 - t241 * t265; t220, t239 * t259 - t253 * t292, -t233 * t259, 0, t224, 0; t223, t237 * t259 - t251 * t292, -t229 * t259, 0, -t222, 0; 0, t248 * t259 - t260 * t280, -t241 * t259, 0, t228, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end