% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:56
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR15_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR15_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
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
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
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
	% StartTime: 2019-10-10 12:56:26
	% EndTime: 2019-10-10 12:56:26
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
	% StartTime: 2019-10-10 12:56:27
	% EndTime: 2019-10-10 12:56:27
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
	% StartTime: 2019-10-10 12:56:27
	% EndTime: 2019-10-10 12:56:27
	% DurationCPUTime: 0.30s
	% Computational Cost: add. (136->42), mult. (408->90), div. (0->0), fcn. (571->12), ass. (0->50)
	t209 = sin(qJ(2));
	t210 = sin(qJ(1));
	t213 = cos(qJ(2));
	t214 = cos(qJ(1));
	t233 = cos(pkin(6));
	t217 = t214 * t233;
	t198 = t209 * t217 + t210 * t213;
	t208 = sin(qJ(3));
	t212 = cos(qJ(3));
	t197 = t210 * t209 - t213 * t217;
	t204 = sin(pkin(7));
	t206 = cos(pkin(7));
	t205 = sin(pkin(6));
	t227 = t205 * t214;
	t215 = t197 * t206 + t204 * t227;
	t184 = -t198 * t212 + t215 * t208;
	t191 = -t197 * t204 + t206 * t227;
	t207 = sin(qJ(4));
	t211 = cos(qJ(4));
	t237 = t184 * t207 - t191 * t211;
	t236 = -t184 * t211 - t191 * t207;
	t231 = t204 * t205;
	t230 = t204 * t207;
	t229 = t204 * t211;
	t228 = t205 * t210;
	t226 = t206 * t208;
	t225 = t206 * t212;
	t224 = t208 * t209;
	t223 = t208 * t213;
	t222 = t209 * t212;
	t221 = t212 * t213;
	t220 = t209 * t231;
	t219 = t204 * t228;
	t218 = t210 * t233;
	t216 = t233 * t204;
	t182 = -t198 * t208 - t215 * t212;
	t200 = -t209 * t218 + t214 * t213;
	t199 = -t214 * t209 - t213 * t218;
	t196 = t233 * t206 - t213 * t231;
	t195 = (-t206 * t224 + t221) * t205;
	t193 = -t199 * t204 + t206 * t228;
	t190 = t208 * t216 + (t206 * t223 + t222) * t205;
	t189 = t212 * t216 + (t206 * t221 - t224) * t205;
	t188 = t199 * t212 - t200 * t226;
	t187 = -t197 * t212 - t198 * t226;
	t186 = t200 * t212 + (t199 * t206 + t219) * t208;
	t185 = -t199 * t225 + t200 * t208 - t212 * t219;
	t181 = t186 * t211 + t193 * t207;
	t180 = t186 * t207 - t193 * t211;
	t1 = [t182, t199 * t208 + t200 * t225, t186, 0, 0, 0; t185, -t197 * t208 + t198 * t225, -t184, 0, 0, 0; 0, (t206 * t222 + t223) * t205, t190, 0, 0, 0; t236, -t188 * t211 - t200 * t230, t185 * t211, t180, 0, 0; -t181, -t187 * t211 - t198 * t230, -t182 * t211, -t237, 0, 0; 0, -t195 * t211 - t207 * t220, -t189 * t211, t190 * t207 - t196 * t211, 0, 0; t237, t188 * t207 - t200 * t229, -t185 * t207, t181, 0, 0; t180, t187 * t207 - t198 * t229, t182 * t207, t236, 0, 0; 0, t195 * t207 - t211 * t220, t189 * t207, t190 * t211 + t196 * t207, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:56:27
	% EndTime: 2019-10-10 12:56:28
	% DurationCPUTime: 0.56s
	% Computational Cost: add. (290->66), mult. (867->135), div. (0->0), fcn. (1196->14), ass. (0->69)
	t252 = sin(qJ(2));
	t253 = sin(qJ(1));
	t257 = cos(qJ(2));
	t258 = cos(qJ(1));
	t281 = cos(pkin(6));
	t261 = t258 * t281;
	t238 = t252 * t261 + t253 * t257;
	t251 = sin(qJ(3));
	t256 = cos(qJ(3));
	t237 = t253 * t252 - t257 * t261;
	t248 = cos(pkin(7));
	t246 = sin(pkin(7));
	t247 = sin(pkin(6));
	t274 = t247 * t258;
	t263 = t246 * t274;
	t259 = t237 * t248 + t263;
	t219 = -t238 * t256 + t259 * t251;
	t229 = -t237 * t246 + t248 * t274;
	t250 = sin(qJ(4));
	t255 = cos(qJ(4));
	t209 = t219 * t250 - t229 * t255;
	t249 = sin(qJ(6));
	t286 = t209 * t249;
	t254 = cos(qJ(6));
	t285 = t209 * t254;
	t284 = t219 * t255 + t229 * t250;
	t280 = t238 * t251;
	t278 = t246 * t247;
	t277 = t246 * t250;
	t276 = t246 * t255;
	t275 = t247 * t253;
	t273 = t248 * t251;
	t272 = t248 * t256;
	t271 = t249 * t250;
	t270 = t250 * t254;
	t269 = t251 * t252;
	t268 = t251 * t257;
	t267 = t252 * t256;
	t266 = t256 * t257;
	t265 = t252 * t278;
	t264 = t246 * t275;
	t262 = t253 * t281;
	t260 = t281 * t246;
	t240 = -t252 * t262 + t258 * t257;
	t239 = -t258 * t252 - t257 * t262;
	t236 = t281 * t248 - t257 * t278;
	t235 = (-t248 * t269 + t266) * t247;
	t234 = (t248 * t267 + t268) * t247;
	t231 = -t239 * t246 + t248 * t275;
	t228 = t251 * t260 + (t248 * t268 + t267) * t247;
	t227 = -t256 * t260 + (-t248 * t266 + t269) * t247;
	t226 = t235 * t250 - t255 * t265;
	t225 = t239 * t256 - t240 * t273;
	t224 = t239 * t251 + t240 * t272;
	t223 = -t237 * t256 - t238 * t273;
	t222 = -t237 * t251 + t238 * t272;
	t221 = t240 * t256 + (t239 * t248 + t264) * t251;
	t220 = -t239 * t272 + t240 * t251 - t256 * t264;
	t218 = -t259 * t256 - t280;
	t216 = t237 * t272 + t256 * t263 + t280;
	t215 = t228 * t255 + t236 * t250;
	t214 = t228 * t250 - t236 * t255;
	t213 = t225 * t250 - t240 * t276;
	t212 = t223 * t250 - t238 * t276;
	t211 = t221 * t255 + t231 * t250;
	t210 = t221 * t250 - t231 * t255;
	t206 = t210 * t249 + t220 * t254;
	t205 = t210 * t254 - t220 * t249;
	t1 = [t218 * t254 + t286, t213 * t249 + t224 * t254, -t220 * t271 + t221 * t254, t211 * t249, 0, t205; t206, t212 * t249 + t222 * t254, -t216 * t271 - t219 * t254, -t284 * t249, 0, -t216 * t249 - t285; 0, t226 * t249 + t234 * t254, -t227 * t271 + t228 * t254, t215 * t249, 0, t214 * t254 - t227 * t249; -t218 * t249 + t285, t213 * t254 - t224 * t249, -t220 * t270 - t221 * t249, t211 * t254, 0, -t206; t205, t212 * t254 - t222 * t249, -t216 * t270 + t219 * t249, -t284 * t254, 0, -t216 * t254 + t286; 0, t226 * t254 - t234 * t249, -t227 * t270 - t228 * t249, t215 * t254, 0, -t214 * t249 - t227 * t254; t284, t225 * t255 + t240 * t277, -t220 * t255, -t210, 0, 0; t211, t223 * t255 + t238 * t277, -t216 * t255, t209, 0, 0; 0, t235 * t255 + t250 * t265, -t227 * t255, -t214, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end