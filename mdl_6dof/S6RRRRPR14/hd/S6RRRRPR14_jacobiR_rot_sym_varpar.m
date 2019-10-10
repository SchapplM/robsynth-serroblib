% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR14
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
% Datum: 2019-10-10 12:54
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR14_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR14_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
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
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
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
	% StartTime: 2019-10-10 12:54:13
	% EndTime: 2019-10-10 12:54:13
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
	% StartTime: 2019-10-10 12:54:13
	% EndTime: 2019-10-10 12:54:13
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
	% StartTime: 2019-10-10 12:54:14
	% EndTime: 2019-10-10 12:54:14
	% DurationCPUTime: 0.45s
	% Computational Cost: add. (231->57), mult. (689->127), div. (0->0), fcn. (950->14), ass. (0->61)
	t235 = sin(qJ(2));
	t236 = sin(qJ(1));
	t239 = cos(qJ(2));
	t240 = cos(qJ(1));
	t263 = cos(pkin(6));
	t245 = t240 * t263;
	t222 = t235 * t245 + t236 * t239;
	t234 = sin(qJ(3));
	t238 = cos(qJ(3));
	t221 = t236 * t235 - t239 * t245;
	t229 = sin(pkin(7));
	t232 = cos(pkin(7));
	t230 = sin(pkin(6));
	t256 = t230 * t240;
	t243 = t221 * t232 + t229 * t256;
	t204 = -t222 * t238 + t243 * t234;
	t215 = -t221 * t229 + t232 * t256;
	t233 = sin(qJ(4));
	t237 = cos(qJ(4));
	t194 = t204 * t233 - t215 * t237;
	t195 = t204 * t237 + t215 * t233;
	t228 = sin(pkin(13));
	t261 = t228 * t237;
	t260 = t229 * t230;
	t259 = t229 * t233;
	t258 = t229 * t237;
	t257 = t230 * t236;
	t231 = cos(pkin(13));
	t255 = t231 * t237;
	t254 = t232 * t234;
	t253 = t232 * t238;
	t252 = t234 * t235;
	t251 = t234 * t239;
	t250 = t235 * t238;
	t249 = t238 * t239;
	t248 = t235 * t260;
	t247 = t229 * t257;
	t246 = t236 * t263;
	t244 = t263 * t229;
	t223 = -t240 * t235 - t239 * t246;
	t242 = -t223 * t229 + t232 * t257;
	t241 = -t222 * t234 - t243 * t238;
	t224 = -t235 * t246 + t240 * t239;
	t220 = t263 * t232 - t239 * t260;
	t219 = (-t232 * t252 + t249) * t230;
	t218 = (t232 * t250 + t251) * t230;
	t214 = t234 * t244 + (t232 * t251 + t250) * t230;
	t213 = t238 * t244 + (t232 * t249 - t252) * t230;
	t211 = t219 * t237 + t233 * t248;
	t210 = t223 * t238 - t224 * t254;
	t209 = t223 * t234 + t224 * t253;
	t208 = -t221 * t238 - t222 * t254;
	t207 = -t221 * t234 + t222 * t253;
	t206 = t224 * t238 + (t223 * t232 + t247) * t234;
	t205 = -t223 * t253 + t224 * t234 - t238 * t247;
	t200 = -t214 * t233 + t220 * t237;
	t199 = t210 * t237 + t224 * t259;
	t198 = t208 * t237 + t222 * t259;
	t197 = t206 * t237 + t242 * t233;
	t196 = t206 * t233 - t242 * t237;
	t1 = [t195 * t231 + t228 * t241, t199 * t231 + t209 * t228, -t205 * t255 + t206 * t228, -t196 * t231, 0, 0; t197 * t231 + t205 * t228, t198 * t231 + t207 * t228, -t204 * t228 + t241 * t255, t194 * t231, 0, 0; 0, t211 * t231 + t218 * t228, t213 * t255 + t214 * t228, t200 * t231, 0, 0; -t195 * t228 + t231 * t241, -t199 * t228 + t209 * t231, t205 * t261 + t206 * t231, t196 * t228, 0, 0; -t197 * t228 + t205 * t231, -t198 * t228 + t207 * t231, -t204 * t231 - t241 * t261, -t194 * t228, 0, 0; 0, -t211 * t228 + t218 * t231, -t213 * t261 + t214 * t231, -t200 * t228, 0, 0; t194, t210 * t233 - t224 * t258, -t205 * t233, t197, 0, 0; t196, t208 * t233 - t222 * t258, t241 * t233, -t195, 0, 0; 0, t219 * t233 - t237 * t248, t213 * t233, t214 * t237 + t220 * t233, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:14
	% EndTime: 2019-10-10 12:54:14
	% DurationCPUTime: 0.52s
	% Computational Cost: add. (343->67), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->70)
	t260 = sin(qJ(2));
	t261 = sin(qJ(1));
	t264 = cos(qJ(2));
	t265 = cos(qJ(1));
	t289 = cos(pkin(6));
	t269 = t265 * t289;
	t244 = t260 * t269 + t261 * t264;
	t259 = sin(qJ(3));
	t263 = cos(qJ(3));
	t243 = t261 * t260 - t264 * t269;
	t257 = cos(pkin(7));
	t255 = sin(pkin(7));
	t256 = sin(pkin(6));
	t280 = t256 * t265;
	t271 = t255 * t280;
	t267 = t243 * t257 + t271;
	t225 = -t244 * t263 + t267 * t259;
	t236 = -t243 * t255 + t257 * t280;
	t258 = sin(qJ(4));
	t262 = cos(qJ(4));
	t215 = t225 * t262 + t236 * t258;
	t254 = pkin(13) + qJ(6);
	t252 = sin(t254);
	t293 = t215 * t252;
	t253 = cos(t254);
	t292 = t215 * t253;
	t213 = t225 * t258 - t236 * t262;
	t288 = t244 * t259;
	t286 = t252 * t262;
	t285 = t253 * t262;
	t284 = t255 * t256;
	t283 = t255 * t258;
	t282 = t255 * t262;
	t281 = t256 * t261;
	t279 = t257 * t259;
	t278 = t257 * t263;
	t277 = t259 * t260;
	t276 = t259 * t264;
	t275 = t260 * t263;
	t274 = t263 * t264;
	t273 = t260 * t284;
	t272 = t255 * t281;
	t270 = t261 * t289;
	t268 = t289 * t255;
	t245 = -t265 * t260 - t264 * t270;
	t266 = -t245 * t255 + t257 * t281;
	t246 = -t260 * t270 + t265 * t264;
	t242 = t289 * t257 - t264 * t284;
	t241 = (-t257 * t277 + t274) * t256;
	t240 = (t257 * t275 + t276) * t256;
	t235 = t259 * t268 + (t257 * t276 + t275) * t256;
	t234 = -t263 * t268 + (-t257 * t274 + t277) * t256;
	t232 = t241 * t262 + t258 * t273;
	t231 = t245 * t263 - t246 * t279;
	t230 = t245 * t259 + t246 * t278;
	t229 = -t243 * t263 - t244 * t279;
	t228 = -t243 * t259 + t244 * t278;
	t227 = t246 * t263 + (t245 * t257 + t272) * t259;
	t226 = -t245 * t278 + t246 * t259 - t263 * t272;
	t224 = -t267 * t263 - t288;
	t222 = t243 * t278 + t263 * t271 + t288;
	t221 = t235 * t262 + t242 * t258;
	t220 = -t235 * t258 + t242 * t262;
	t219 = t231 * t262 + t246 * t283;
	t218 = t229 * t262 + t244 * t283;
	t217 = t227 * t262 + t266 * t258;
	t216 = t227 * t258 - t266 * t262;
	t212 = t217 * t253 + t226 * t252;
	t211 = -t217 * t252 + t226 * t253;
	t1 = [t224 * t252 + t292, t219 * t253 + t230 * t252, -t226 * t285 + t227 * t252, -t216 * t253, 0, t211; t212, t218 * t253 + t228 * t252, -t222 * t285 - t225 * t252, t213 * t253, 0, t222 * t253 + t293; 0, t232 * t253 + t240 * t252, -t234 * t285 + t235 * t252, t220 * t253, 0, -t221 * t252 + t234 * t253; t224 * t253 - t293, -t219 * t252 + t230 * t253, t226 * t286 + t227 * t253, t216 * t252, 0, -t212; t211, -t218 * t252 + t228 * t253, t222 * t286 - t225 * t253, -t213 * t252, 0, -t222 * t252 + t292; 0, -t232 * t252 + t240 * t253, t234 * t286 + t235 * t253, -t220 * t252, 0, -t221 * t253 - t234 * t252; t213, t231 * t258 - t246 * t282, -t226 * t258, t217, 0, 0; t216, t229 * t258 - t244 * t282, -t222 * t258, -t215, 0, 0; 0, t241 * t258 - t262 * t273, -t234 * t258, t221, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end