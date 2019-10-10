% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR12_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR12_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
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
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:58
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
	% StartTime: 2019-10-10 12:49:58
	% EndTime: 2019-10-10 12:49:59
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
	% StartTime: 2019-10-10 12:49:59
	% EndTime: 2019-10-10 12:49:59
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
	% StartTime: 2019-10-10 12:49:59
	% EndTime: 2019-10-10 12:49:59
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (174->43), mult. (408->90), div. (0->0), fcn. (571->12), ass. (0->51)
	t194 = sin(qJ(2));
	t195 = sin(qJ(1));
	t197 = cos(qJ(2));
	t198 = cos(qJ(1));
	t217 = cos(pkin(6));
	t201 = t198 * t217;
	t181 = t194 * t201 + t195 * t197;
	t193 = sin(qJ(3));
	t196 = cos(qJ(3));
	t180 = t195 * t194 - t197 * t201;
	t190 = sin(pkin(7));
	t192 = cos(pkin(7));
	t191 = sin(pkin(6));
	t211 = t191 * t198;
	t199 = t180 * t192 + t190 * t211;
	t167 = -t181 * t196 + t199 * t193;
	t174 = -t180 * t190 + t192 * t211;
	t189 = qJ(4) + pkin(13);
	t187 = sin(t189);
	t188 = cos(t189);
	t221 = t167 * t187 - t174 * t188;
	t220 = t167 * t188 + t174 * t187;
	t215 = t187 * t190;
	t214 = t188 * t190;
	t213 = t190 * t191;
	t212 = t191 * t195;
	t210 = t192 * t193;
	t209 = t192 * t196;
	t208 = t193 * t194;
	t207 = t193 * t197;
	t206 = t194 * t196;
	t205 = t196 * t197;
	t204 = t194 * t213;
	t203 = t190 * t212;
	t202 = t195 * t217;
	t200 = t217 * t190;
	t165 = -t181 * t193 - t199 * t196;
	t183 = -t194 * t202 + t198 * t197;
	t182 = -t198 * t194 - t197 * t202;
	t179 = t217 * t192 - t197 * t213;
	t178 = (-t192 * t208 + t205) * t191;
	t176 = -t182 * t190 + t192 * t212;
	t173 = t193 * t200 + (t192 * t207 + t206) * t191;
	t172 = t196 * t200 + (t192 * t205 - t208) * t191;
	t171 = t182 * t196 - t183 * t210;
	t170 = -t180 * t196 - t181 * t210;
	t169 = t183 * t196 + (t182 * t192 + t203) * t193;
	t168 = -t182 * t209 + t183 * t193 - t196 * t203;
	t164 = t169 * t188 + t176 * t187;
	t163 = -t169 * t187 + t176 * t188;
	t1 = [t220, t171 * t188 + t183 * t215, -t168 * t188, t163, 0, 0; t164, t170 * t188 + t181 * t215, t165 * t188, t221, 0, 0; 0, t178 * t188 + t187 * t204, t172 * t188, -t173 * t187 + t179 * t188, 0, 0; -t221, -t171 * t187 + t183 * t214, t168 * t187, -t164, 0, 0; t163, -t170 * t187 + t181 * t214, -t165 * t187, t220, 0, 0; 0, -t178 * t187 + t188 * t204, -t172 * t187, -t173 * t188 - t179 * t187, 0, 0; t165, t182 * t193 + t183 * t209, t169, 0, 0, 0; t168, -t180 * t193 + t181 * t209, -t167, 0, 0, 0; 0, (t192 * t206 + t207) * t191, t173, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:50:00
	% EndTime: 2019-10-10 12:50:00
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (362->67), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->70)
	t268 = sin(qJ(2));
	t269 = sin(qJ(1));
	t272 = cos(qJ(2));
	t273 = cos(qJ(1));
	t297 = cos(pkin(6));
	t277 = t273 * t297;
	t252 = t268 * t277 + t269 * t272;
	t267 = sin(qJ(3));
	t271 = cos(qJ(3));
	t251 = t269 * t268 - t272 * t277;
	t265 = cos(pkin(7));
	t263 = sin(pkin(7));
	t264 = sin(pkin(6));
	t288 = t264 * t273;
	t279 = t263 * t288;
	t275 = t251 * t265 + t279;
	t233 = -t252 * t271 + t275 * t267;
	t244 = -t251 * t263 + t265 * t288;
	t262 = qJ(4) + pkin(13);
	t260 = sin(t262);
	t261 = cos(t262);
	t223 = t233 * t261 + t244 * t260;
	t266 = sin(qJ(6));
	t301 = t223 * t266;
	t270 = cos(qJ(6));
	t300 = t223 * t270;
	t221 = t233 * t260 - t244 * t261;
	t296 = t252 * t267;
	t294 = t260 * t263;
	t293 = t261 * t263;
	t292 = t261 * t266;
	t291 = t261 * t270;
	t290 = t263 * t264;
	t289 = t264 * t269;
	t287 = t265 * t267;
	t286 = t265 * t271;
	t285 = t267 * t268;
	t284 = t267 * t272;
	t283 = t268 * t271;
	t282 = t271 * t272;
	t281 = t268 * t290;
	t280 = t263 * t289;
	t278 = t269 * t297;
	t276 = t297 * t263;
	t253 = -t273 * t268 - t272 * t278;
	t274 = -t253 * t263 + t265 * t289;
	t254 = -t268 * t278 + t273 * t272;
	t250 = t297 * t265 - t272 * t290;
	t249 = (-t265 * t285 + t282) * t264;
	t248 = (t265 * t283 + t284) * t264;
	t243 = t267 * t276 + (t265 * t284 + t283) * t264;
	t242 = -t271 * t276 + (-t265 * t282 + t285) * t264;
	t240 = t253 * t271 - t254 * t287;
	t239 = t253 * t267 + t254 * t286;
	t238 = -t251 * t271 - t252 * t287;
	t237 = -t251 * t267 + t252 * t286;
	t236 = t249 * t261 + t260 * t281;
	t235 = t254 * t271 + (t253 * t265 + t280) * t267;
	t234 = -t253 * t286 + t254 * t267 - t271 * t280;
	t232 = -t275 * t271 - t296;
	t230 = t251 * t286 + t271 * t279 + t296;
	t229 = t243 * t261 + t250 * t260;
	t228 = -t243 * t260 + t250 * t261;
	t227 = t240 * t261 + t254 * t294;
	t226 = t238 * t261 + t252 * t294;
	t225 = t235 * t261 + t274 * t260;
	t224 = t235 * t260 - t274 * t261;
	t220 = t225 * t270 + t234 * t266;
	t219 = -t225 * t266 + t234 * t270;
	t1 = [t232 * t266 + t300, t227 * t270 + t239 * t266, -t234 * t291 + t235 * t266, -t224 * t270, 0, t219; t220, t226 * t270 + t237 * t266, -t230 * t291 - t233 * t266, t221 * t270, 0, t230 * t270 + t301; 0, t236 * t270 + t248 * t266, -t242 * t291 + t243 * t266, t228 * t270, 0, -t229 * t266 + t242 * t270; t232 * t270 - t301, -t227 * t266 + t239 * t270, t234 * t292 + t235 * t270, t224 * t266, 0, -t220; t219, -t226 * t266 + t237 * t270, t230 * t292 - t233 * t270, -t221 * t266, 0, -t230 * t266 + t300; 0, -t236 * t266 + t248 * t270, t242 * t292 + t243 * t270, -t228 * t266, 0, -t229 * t270 - t242 * t266; t221, t240 * t260 - t254 * t293, -t234 * t260, t225, 0, 0; t224, t238 * t260 - t252 * t293, -t230 * t260, -t223, 0, 0; 0, t249 * t260 - t261 * t281, -t242 * t260, t229, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end