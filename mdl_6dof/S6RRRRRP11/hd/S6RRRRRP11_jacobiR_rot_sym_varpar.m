% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRRP11
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:13
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP11_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP11_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:13:34
	% EndTime: 2019-10-10 13:13:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:13:34
	% EndTime: 2019-10-10 13:13:34
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
	% StartTime: 2019-10-10 13:13:34
	% EndTime: 2019-10-10 13:13:34
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
	% StartTime: 2019-10-10 13:13:34
	% EndTime: 2019-10-10 13:13:35
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
	% StartTime: 2019-10-10 13:13:35
	% EndTime: 2019-10-10 13:13:35
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
	% StartTime: 2019-10-10 13:13:36
	% EndTime: 2019-10-10 13:13:36
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (293->66), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->69)
	t256 = sin(qJ(2));
	t257 = sin(qJ(1));
	t261 = cos(qJ(2));
	t262 = cos(qJ(1));
	t286 = cos(pkin(6));
	t266 = t262 * t286;
	t242 = t256 * t266 + t257 * t261;
	t255 = sin(qJ(3));
	t260 = cos(qJ(3));
	t241 = t257 * t256 - t261 * t266;
	t252 = cos(pkin(7));
	t250 = sin(pkin(7));
	t251 = sin(pkin(6));
	t279 = t251 * t262;
	t268 = t250 * t279;
	t264 = t241 * t252 + t268;
	t223 = -t242 * t260 + t264 * t255;
	t234 = -t241 * t250 + t252 * t279;
	t254 = sin(qJ(4));
	t259 = cos(qJ(4));
	t213 = t223 * t259 + t234 * t254;
	t253 = sin(qJ(5));
	t290 = t213 * t253;
	t258 = cos(qJ(5));
	t289 = t213 * t258;
	t211 = t223 * t254 - t234 * t259;
	t285 = t242 * t255;
	t283 = t250 * t251;
	t282 = t250 * t254;
	t281 = t250 * t259;
	t280 = t251 * t257;
	t278 = t252 * t255;
	t277 = t252 * t260;
	t276 = t253 * t259;
	t275 = t255 * t256;
	t274 = t255 * t261;
	t273 = t256 * t260;
	t272 = t258 * t259;
	t271 = t260 * t261;
	t270 = t256 * t283;
	t269 = t250 * t280;
	t267 = t257 * t286;
	t265 = t286 * t250;
	t243 = -t262 * t256 - t261 * t267;
	t263 = -t243 * t250 + t252 * t280;
	t244 = -t256 * t267 + t262 * t261;
	t240 = t286 * t252 - t261 * t283;
	t239 = (-t252 * t275 + t271) * t251;
	t238 = (t252 * t273 + t274) * t251;
	t233 = t255 * t265 + (t252 * t274 + t273) * t251;
	t232 = -t260 * t265 + (-t252 * t271 + t275) * t251;
	t230 = t239 * t259 + t254 * t270;
	t229 = t243 * t260 - t244 * t278;
	t228 = t243 * t255 + t244 * t277;
	t227 = -t241 * t260 - t242 * t278;
	t226 = -t241 * t255 + t242 * t277;
	t225 = t244 * t260 + (t243 * t252 + t269) * t255;
	t224 = -t243 * t277 + t244 * t255 - t260 * t269;
	t222 = -t264 * t260 - t285;
	t220 = t241 * t277 + t260 * t268 + t285;
	t219 = t233 * t259 + t240 * t254;
	t218 = -t233 * t254 + t240 * t259;
	t217 = t229 * t259 + t244 * t282;
	t216 = t227 * t259 + t242 * t282;
	t215 = t225 * t259 + t263 * t254;
	t214 = t225 * t254 - t263 * t259;
	t210 = t215 * t258 + t224 * t253;
	t209 = -t215 * t253 + t224 * t258;
	t1 = [t222 * t253 + t289, t217 * t258 + t228 * t253, -t224 * t272 + t225 * t253, -t214 * t258, t209, 0; t210, t216 * t258 + t226 * t253, -t220 * t272 - t223 * t253, t211 * t258, t220 * t258 + t290, 0; 0, t230 * t258 + t238 * t253, -t232 * t272 + t233 * t253, t218 * t258, -t219 * t253 + t232 * t258, 0; t222 * t258 - t290, -t217 * t253 + t228 * t258, t224 * t276 + t225 * t258, t214 * t253, -t210, 0; t209, -t216 * t253 + t226 * t258, t220 * t276 - t223 * t258, -t211 * t253, -t220 * t253 + t289, 0; 0, -t230 * t253 + t238 * t258, t232 * t276 + t233 * t258, -t218 * t253, -t219 * t258 - t232 * t253, 0; t211, t229 * t254 - t244 * t281, -t224 * t254, t215, 0, 0; t214, t227 * t254 - t242 * t281, -t220 * t254, -t213, 0, 0; 0, t239 * t254 - t259 * t270, -t232 * t254, t219, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:13:35
	% EndTime: 2019-10-10 13:13:36
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (293->66), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->69)
	t257 = sin(qJ(2));
	t258 = sin(qJ(1));
	t262 = cos(qJ(2));
	t263 = cos(qJ(1));
	t287 = cos(pkin(6));
	t267 = t263 * t287;
	t243 = t257 * t267 + t258 * t262;
	t256 = sin(qJ(3));
	t261 = cos(qJ(3));
	t242 = t258 * t257 - t262 * t267;
	t253 = cos(pkin(7));
	t251 = sin(pkin(7));
	t252 = sin(pkin(6));
	t280 = t252 * t263;
	t269 = t251 * t280;
	t265 = t242 * t253 + t269;
	t224 = -t243 * t261 + t265 * t256;
	t235 = -t242 * t251 + t253 * t280;
	t255 = sin(qJ(4));
	t260 = cos(qJ(4));
	t214 = t224 * t260 + t235 * t255;
	t254 = sin(qJ(5));
	t291 = t214 * t254;
	t259 = cos(qJ(5));
	t290 = t214 * t259;
	t212 = t224 * t255 - t235 * t260;
	t286 = t243 * t256;
	t284 = t251 * t252;
	t283 = t251 * t255;
	t282 = t251 * t260;
	t281 = t252 * t258;
	t279 = t253 * t256;
	t278 = t253 * t261;
	t277 = t254 * t260;
	t276 = t256 * t257;
	t275 = t256 * t262;
	t274 = t257 * t261;
	t273 = t259 * t260;
	t272 = t261 * t262;
	t271 = t257 * t284;
	t270 = t251 * t281;
	t268 = t258 * t287;
	t266 = t287 * t251;
	t244 = -t263 * t257 - t262 * t268;
	t264 = -t244 * t251 + t253 * t281;
	t245 = -t257 * t268 + t263 * t262;
	t241 = t287 * t253 - t262 * t284;
	t240 = (-t253 * t276 + t272) * t252;
	t239 = (t253 * t274 + t275) * t252;
	t234 = t256 * t266 + (t253 * t275 + t274) * t252;
	t233 = -t261 * t266 + (-t253 * t272 + t276) * t252;
	t231 = t240 * t260 + t255 * t271;
	t230 = t244 * t261 - t245 * t279;
	t229 = t244 * t256 + t245 * t278;
	t228 = -t242 * t261 - t243 * t279;
	t227 = -t242 * t256 + t243 * t278;
	t226 = t245 * t261 + (t244 * t253 + t270) * t256;
	t225 = -t244 * t278 + t245 * t256 - t261 * t270;
	t223 = -t265 * t261 - t286;
	t221 = t242 * t278 + t261 * t269 + t286;
	t220 = t234 * t260 + t241 * t255;
	t219 = -t234 * t255 + t241 * t260;
	t218 = t230 * t260 + t245 * t283;
	t217 = t228 * t260 + t243 * t283;
	t216 = t226 * t260 + t264 * t255;
	t215 = t226 * t255 - t264 * t260;
	t211 = t216 * t259 + t225 * t254;
	t210 = -t216 * t254 + t225 * t259;
	t1 = [t223 * t254 + t290, t218 * t259 + t229 * t254, -t225 * t273 + t226 * t254, -t215 * t259, t210, 0; t211, t217 * t259 + t227 * t254, -t221 * t273 - t224 * t254, t212 * t259, t221 * t259 + t291, 0; 0, t231 * t259 + t239 * t254, -t233 * t273 + t234 * t254, t219 * t259, -t220 * t254 + t233 * t259, 0; t223 * t259 - t291, -t218 * t254 + t229 * t259, t225 * t277 + t226 * t259, t215 * t254, -t211, 0; t210, -t217 * t254 + t227 * t259, t221 * t277 - t224 * t259, -t212 * t254, -t221 * t254 + t290, 0; 0, -t231 * t254 + t239 * t259, t233 * t277 + t234 * t259, -t219 * t254, -t220 * t259 - t233 * t254, 0; t212, t230 * t255 - t245 * t282, -t225 * t255, t216, 0, 0; t215, t228 * t255 - t243 * t282, -t221 * t255, -t214, 0, 0; 0, t240 * t255 - t260 * t271, -t233 * t255, t220, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end