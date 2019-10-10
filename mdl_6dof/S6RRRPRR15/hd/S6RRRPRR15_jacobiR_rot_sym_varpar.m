% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR15_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR15_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
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
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
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
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
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
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (58->24), mult. (178->57), div. (0->0), fcn. (255->10), ass. (0->32)
	t143 = sin(qJ(2));
	t144 = sin(qJ(1));
	t146 = cos(qJ(2));
	t147 = cos(qJ(1));
	t163 = cos(pkin(6));
	t152 = t147 * t163;
	t134 = t143 * t152 + t144 * t146;
	t142 = sin(qJ(3));
	t145 = cos(qJ(3));
	t133 = t144 * t143 - t146 * t152;
	t139 = sin(pkin(7));
	t141 = cos(pkin(7));
	t140 = sin(pkin(6));
	t160 = t140 * t147;
	t150 = t133 * t141 + t139 * t160;
	t164 = t134 * t142 + t150 * t145;
	t161 = t140 * t144;
	t159 = t141 * t142;
	t158 = t141 * t145;
	t157 = t142 * t143;
	t156 = t142 * t146;
	t155 = t143 * t145;
	t154 = t145 * t146;
	t153 = t144 * t163;
	t151 = t163 * t139;
	t135 = -t147 * t143 - t146 * t153;
	t149 = t135 * t141 + t139 * t161;
	t148 = t134 * t145 - t150 * t142;
	t136 = -t143 * t153 + t147 * t146;
	t132 = t136 * t145 + t149 * t142;
	t131 = t136 * t142 - t149 * t145;
	t1 = [-t133 * t139 + t141 * t160, t136 * t139, 0, 0, 0, 0; -t135 * t139 + t141 * t161, t134 * t139, 0, 0, 0, 0; 0, t140 * t143 * t139, 0, 0, 0, 0; t148, -t135 * t145 + t136 * t159, t131, 0, 0, 0; -t132, t133 * t145 + t134 * t159, t164, 0, 0, 0; 0, (t141 * t157 - t154) * t140, -t145 * t151 + (-t141 * t154 + t157) * t140, 0, 0, 0; -t164, t135 * t142 + t136 * t158, t132, 0, 0, 0; t131, -t133 * t142 + t134 * t158, t148, 0, 0, 0; 0, (t141 * t155 + t156) * t140, t142 * t151 + (t141 * t156 + t155) * t140, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:24
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (133->42), mult. (408->89), div. (0->0), fcn. (571->12), ass. (0->50)
	t182 = sin(qJ(2));
	t183 = sin(qJ(1));
	t186 = cos(qJ(2));
	t187 = cos(qJ(1));
	t207 = cos(pkin(6));
	t191 = t187 * t207;
	t172 = t182 * t191 + t183 * t186;
	t181 = sin(qJ(3));
	t185 = cos(qJ(3));
	t171 = t183 * t182 - t186 * t191;
	t177 = sin(pkin(7));
	t179 = cos(pkin(7));
	t178 = sin(pkin(6));
	t200 = t178 * t187;
	t189 = t171 * t179 + t177 * t200;
	t157 = t172 * t181 + t189 * t185;
	t166 = -t171 * t177 + t179 * t200;
	t180 = sin(qJ(5));
	t184 = cos(qJ(5));
	t212 = -t157 * t180 + t166 * t184;
	t211 = t157 * t184 + t166 * t180;
	t208 = -t172 * t185 + t189 * t181;
	t204 = t177 * t178;
	t203 = t177 * t180;
	t202 = t177 * t184;
	t201 = t178 * t183;
	t199 = t179 * t181;
	t198 = t179 * t185;
	t197 = t181 * t182;
	t196 = t181 * t186;
	t195 = t182 * t185;
	t194 = t185 * t186;
	t193 = t182 * t204;
	t192 = t183 * t207;
	t190 = t207 * t177;
	t173 = -t187 * t182 - t186 * t192;
	t188 = t173 * t179 + t177 * t201;
	t174 = -t182 * t192 + t187 * t186;
	t170 = t207 * t179 - t186 * t204;
	t169 = (t179 * t195 + t196) * t178;
	t168 = -t173 * t177 + t179 * t201;
	t165 = t181 * t190 + (t179 * t196 + t195) * t178;
	t164 = -t185 * t190 + (-t179 * t194 + t197) * t178;
	t163 = t173 * t181 + t174 * t198;
	t162 = -t171 * t181 + t172 * t198;
	t161 = t174 * t185 + t188 * t181;
	t160 = t174 * t181 - t188 * t185;
	t156 = t160 * t180 + t168 * t184;
	t155 = t160 * t184 - t168 * t180;
	t1 = [t212, t163 * t180 + t174 * t202, t161 * t180, 0, t155, 0; t156, t162 * t180 + t172 * t202, -t208 * t180, 0, t211, 0; 0, t169 * t180 + t184 * t193, t165 * t180, 0, t164 * t184 - t170 * t180, 0; -t211, t163 * t184 - t174 * t203, t161 * t184, 0, -t156, 0; t155, t162 * t184 - t172 * t203, -t208 * t184, 0, t212, 0; 0, t169 * t184 - t180 * t193, t165 * t184, 0, -t164 * t180 - t170 * t184, 0; t208, t173 * t185 - t174 * t199, -t160, 0, 0, 0; t161, -t171 * t185 - t172 * t199, -t157, 0, 0, 0; 0, (-t179 * t197 + t194) * t178, -t164, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:24
	% EndTime: 2019-10-10 12:18:25
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (296->66), mult. (867->137), div. (0->0), fcn. (1196->14), ass. (0->72)
	t258 = sin(qJ(2));
	t259 = sin(qJ(1));
	t263 = cos(qJ(2));
	t264 = cos(qJ(1));
	t287 = cos(pkin(6));
	t267 = t264 * t287;
	t244 = t258 * t267 + t259 * t263;
	t257 = sin(qJ(3));
	t262 = cos(qJ(3));
	t243 = t259 * t258 - t263 * t267;
	t254 = cos(pkin(7));
	t252 = sin(pkin(7));
	t253 = sin(pkin(6));
	t280 = t253 * t264;
	t269 = t252 * t280;
	t265 = t243 * t254 + t269;
	t225 = -t244 * t262 + t265 * t257;
	t255 = sin(qJ(6));
	t291 = t225 * t255;
	t260 = cos(qJ(6));
	t290 = t225 * t260;
	t235 = -t243 * t252 + t254 * t280;
	t256 = sin(qJ(5));
	t289 = t235 * t256;
	t261 = cos(qJ(5));
	t288 = t235 * t261;
	t286 = t244 * t257;
	t284 = t252 * t253;
	t283 = t252 * t256;
	t282 = t252 * t261;
	t281 = t253 * t259;
	t279 = t254 * t257;
	t278 = t254 * t262;
	t277 = t255 * t256;
	t276 = t256 * t260;
	t275 = t257 * t258;
	t274 = t257 * t263;
	t273 = t258 * t262;
	t272 = t262 * t263;
	t271 = t258 * t284;
	t270 = t252 * t281;
	t268 = t259 * t287;
	t266 = t287 * t252;
	t246 = -t258 * t268 + t264 * t263;
	t245 = -t264 * t258 - t263 * t268;
	t242 = t287 * t254 - t263 * t284;
	t241 = (-t254 * t275 + t272) * t253;
	t240 = (t254 * t273 + t274) * t253;
	t237 = -t245 * t252 + t254 * t281;
	t234 = t257 * t266 + (t254 * t274 + t273) * t253;
	t233 = -t262 * t266 + (-t254 * t272 + t275) * t253;
	t232 = t240 * t256 + t261 * t271;
	t231 = t245 * t262 - t246 * t279;
	t230 = t245 * t257 + t246 * t278;
	t229 = -t243 * t262 - t244 * t279;
	t228 = -t243 * t257 + t244 * t278;
	t227 = t246 * t262 + (t245 * t254 + t270) * t257;
	t226 = -t245 * t278 + t246 * t257 - t262 * t270;
	t224 = -t265 * t262 - t286;
	t222 = t243 * t278 + t262 * t269 + t286;
	t221 = t233 * t256 + t242 * t261;
	t220 = t233 * t261 - t242 * t256;
	t218 = t230 * t256 + t246 * t282;
	t217 = t228 * t256 + t244 * t282;
	t216 = t226 * t256 + t237 * t261;
	t215 = -t226 * t261 + t237 * t256;
	t214 = t222 * t256 - t288;
	t213 = t222 * t261 + t289;
	t212 = t224 * t256 + t288;
	t211 = t216 * t260 + t227 * t255;
	t210 = -t216 * t255 + t227 * t260;
	t1 = [t212 * t260 + t291, t218 * t260 + t231 * t255, -t226 * t255 + t227 * t276, 0, -t215 * t260, t210; t211, t217 * t260 + t229 * t255, -t222 * t255 - t225 * t276, 0, t213 * t260, -t214 * t255 - t290; 0, t232 * t260 + t241 * t255, -t233 * t255 + t234 * t276, 0, t220 * t260, -t221 * t255 + t234 * t260; -t212 * t255 + t290, -t218 * t255 + t231 * t260, -t226 * t260 - t227 * t277, 0, t215 * t255, -t211; t210, -t217 * t255 + t229 * t260, -t222 * t260 + t225 * t277, 0, -t213 * t255, -t214 * t260 + t291; 0, -t232 * t255 + t241 * t260, -t233 * t260 - t234 * t277, 0, -t220 * t255, -t221 * t260 - t234 * t255; -t224 * t261 + t289, -t230 * t261 + t246 * t283, -t227 * t261, 0, t216, 0; t215, -t228 * t261 + t244 * t283, t225 * t261, 0, t214, 0; 0, -t240 * t261 + t256 * t271, -t234 * t261, 0, t221, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end