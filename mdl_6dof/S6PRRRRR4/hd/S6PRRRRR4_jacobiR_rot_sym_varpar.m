% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:19
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR4_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR4_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(13));
	t20 = sin(pkin(6));
	t19 = sin(pkin(13));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:23
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (38->20), mult. (118->51), div. (0->0), fcn. (169->10), ass. (0->27)
	t81 = sin(pkin(7));
	t82 = sin(pkin(6));
	t101 = t81 * t82;
	t85 = cos(pkin(6));
	t100 = t81 * t85;
	t84 = cos(pkin(7));
	t86 = sin(qJ(3));
	t99 = t84 * t86;
	t88 = cos(qJ(3));
	t98 = t84 * t88;
	t87 = sin(qJ(2));
	t97 = t85 * t87;
	t89 = cos(qJ(2));
	t96 = t85 * t89;
	t95 = t86 * t87;
	t94 = t86 * t89;
	t93 = t87 * t88;
	t92 = t88 * t89;
	t80 = sin(pkin(13));
	t83 = cos(pkin(13));
	t75 = -t80 * t87 + t83 * t96;
	t91 = t83 * t101 - t75 * t84;
	t77 = -t80 * t96 - t83 * t87;
	t90 = t80 * t101 + t77 * t84;
	t78 = -t80 * t97 + t83 * t89;
	t76 = t80 * t89 + t83 * t97;
	t1 = [0, t77 * t88 - t78 * t99, -t78 * t86 + t90 * t88, 0, 0, 0; 0, t75 * t88 - t76 * t99, -t76 * t86 - t91 * t88, 0, 0, 0; 0, (-t84 * t95 + t92) * t82, t88 * t100 + (t84 * t92 - t95) * t82, 0, 0, 0; 0, -t77 * t86 - t78 * t98, -t78 * t88 - t90 * t86, 0, 0, 0; 0, -t75 * t86 - t76 * t98, -t76 * t88 + t91 * t86, 0, 0, 0; 0, (-t84 * t93 - t94) * t82, -t86 * t100 + (-t84 * t94 - t93) * t82, 0, 0, 0; 0, t78 * t81, 0, 0, 0, 0; 0, t76 * t81, 0, 0, 0, 0; 0, t87 * t101, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:23
	% EndTime: 2019-10-09 23:19:24
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (100->38), mult. (304->88), div. (0->0), fcn. (425->12), ass. (0->45)
	t152 = sin(pkin(7));
	t153 = sin(pkin(6));
	t178 = t152 * t153;
	t156 = cos(pkin(6));
	t177 = t152 * t156;
	t157 = sin(qJ(4));
	t176 = t152 * t157;
	t160 = cos(qJ(4));
	t175 = t152 * t160;
	t155 = cos(pkin(7));
	t174 = t153 * t155;
	t158 = sin(qJ(3));
	t173 = t155 * t158;
	t161 = cos(qJ(3));
	t172 = t155 * t161;
	t159 = sin(qJ(2));
	t171 = t156 * t159;
	t162 = cos(qJ(2));
	t170 = t156 * t162;
	t169 = t158 * t159;
	t168 = t158 * t162;
	t167 = t159 * t161;
	t166 = t161 * t162;
	t165 = t159 * t178;
	t151 = sin(pkin(13));
	t154 = cos(pkin(13));
	t146 = -t151 * t159 + t154 * t170;
	t164 = t146 * t155 - t154 * t178;
	t148 = -t151 * t170 - t154 * t159;
	t163 = t148 * t155 + t151 * t178;
	t149 = -t151 * t171 + t154 * t162;
	t147 = t151 * t162 + t154 * t171;
	t145 = t156 * t155 - t162 * t178;
	t144 = (-t155 * t169 + t166) * t153;
	t143 = -t148 * t152 + t151 * t174;
	t142 = -t146 * t152 - t154 * t174;
	t141 = t158 * t177 + (t155 * t168 + t167) * t153;
	t140 = t161 * t177 + (t155 * t166 - t169) * t153;
	t139 = t148 * t161 - t149 * t173;
	t138 = t146 * t161 - t147 * t173;
	t137 = t149 * t161 + t163 * t158;
	t136 = -t149 * t158 + t163 * t161;
	t135 = t147 * t161 + t164 * t158;
	t134 = -t147 * t158 + t164 * t161;
	t1 = [0, t139 * t160 + t149 * t176, t136 * t160, -t137 * t157 + t143 * t160, 0, 0; 0, t138 * t160 + t147 * t176, t134 * t160, -t135 * t157 + t142 * t160, 0, 0; 0, t144 * t160 + t157 * t165, t140 * t160, -t141 * t157 + t145 * t160, 0, 0; 0, -t139 * t157 + t149 * t175, -t136 * t157, -t137 * t160 - t143 * t157, 0, 0; 0, -t138 * t157 + t147 * t175, -t134 * t157, -t135 * t160 - t142 * t157, 0, 0; 0, -t144 * t157 + t160 * t165, -t140 * t157, -t141 * t160 - t145 * t157, 0, 0; 0, t148 * t158 + t149 * t172, t137, 0, 0, 0; 0, t146 * t158 + t147 * t172, t135, 0, 0, 0; 0, (t155 * t167 + t168) * t153, t141, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:24
	% EndTime: 2019-10-09 23:19:24
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (178->39), mult. (408->88), div. (0->0), fcn. (571->12), ass. (0->52)
	t176 = qJ(4) + qJ(5);
	t174 = sin(t176);
	t178 = sin(pkin(7));
	t202 = t174 * t178;
	t175 = cos(t176);
	t201 = t175 * t178;
	t179 = sin(pkin(6));
	t200 = t178 * t179;
	t182 = cos(pkin(6));
	t199 = t178 * t182;
	t181 = cos(pkin(7));
	t198 = t179 * t181;
	t183 = sin(qJ(3));
	t197 = t181 * t183;
	t185 = cos(qJ(3));
	t196 = t181 * t185;
	t184 = sin(qJ(2));
	t195 = t182 * t184;
	t186 = cos(qJ(2));
	t194 = t182 * t186;
	t193 = t183 * t184;
	t192 = t183 * t186;
	t191 = t184 * t185;
	t190 = t185 * t186;
	t189 = t184 * t200;
	t177 = sin(pkin(13));
	t180 = cos(pkin(13));
	t169 = -t177 * t184 + t180 * t194;
	t188 = t169 * t181 - t180 * t200;
	t171 = -t177 * t194 - t180 * t184;
	t187 = t171 * t181 + t177 * t200;
	t172 = -t177 * t195 + t180 * t186;
	t170 = t177 * t186 + t180 * t195;
	t168 = t182 * t181 - t186 * t200;
	t167 = (-t181 * t193 + t190) * t179;
	t166 = -t171 * t178 + t177 * t198;
	t165 = -t169 * t178 - t180 * t198;
	t164 = t183 * t199 + (t181 * t192 + t191) * t179;
	t163 = t185 * t199 + (t181 * t190 - t193) * t179;
	t162 = t171 * t185 - t172 * t197;
	t161 = t169 * t185 - t170 * t197;
	t160 = t172 * t185 + t187 * t183;
	t159 = -t172 * t183 + t187 * t185;
	t158 = t170 * t185 + t188 * t183;
	t157 = -t170 * t183 + t188 * t185;
	t156 = -t164 * t175 - t168 * t174;
	t155 = -t164 * t174 + t168 * t175;
	t154 = -t160 * t175 - t166 * t174;
	t153 = -t160 * t174 + t166 * t175;
	t152 = -t158 * t175 - t165 * t174;
	t151 = -t158 * t174 + t165 * t175;
	t1 = [0, t162 * t175 + t172 * t202, t159 * t175, t153, t153, 0; 0, t161 * t175 + t170 * t202, t157 * t175, t151, t151, 0; 0, t167 * t175 + t174 * t189, t163 * t175, t155, t155, 0; 0, -t162 * t174 + t172 * t201, -t159 * t174, t154, t154, 0; 0, -t161 * t174 + t170 * t201, -t157 * t174, t152, t152, 0; 0, -t167 * t174 + t175 * t189, -t163 * t174, t156, t156, 0; 0, t171 * t183 + t172 * t196, t160, 0, 0, 0; 0, t169 * t183 + t170 * t196, t158, 0, 0, 0; 0, (t181 * t191 + t192) * t179, t164, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:19:24
	% EndTime: 2019-10-09 23:19:24
	% DurationCPUTime: 0.28s
	% Computational Cost: add. (363->65), mult. (853->133), div. (0->0), fcn. (1177->14), ass. (0->67)
	t249 = sin(pkin(13));
	t252 = cos(pkin(13));
	t257 = sin(qJ(2));
	t254 = cos(pkin(6));
	t260 = cos(qJ(2));
	t267 = t254 * t260;
	t238 = -t249 * t257 + t252 * t267;
	t268 = t254 * t257;
	t239 = t249 * t260 + t252 * t268;
	t253 = cos(pkin(7));
	t256 = sin(qJ(3));
	t259 = cos(qJ(3));
	t250 = sin(pkin(7));
	t251 = sin(pkin(6));
	t273 = t250 * t251;
	t221 = t239 * t259 + (t238 * t253 - t252 * t273) * t256;
	t271 = t251 * t253;
	t231 = -t238 * t250 - t252 * t271;
	t248 = qJ(4) + qJ(5);
	t246 = sin(t248);
	t247 = cos(t248);
	t211 = -t221 * t246 + t231 * t247;
	t255 = sin(qJ(6));
	t280 = t211 * t255;
	t240 = -t249 * t267 - t252 * t257;
	t241 = -t249 * t268 + t252 * t260;
	t223 = t241 * t259 + (t240 * t253 + t249 * t273) * t256;
	t232 = -t240 * t250 + t249 * t271;
	t213 = -t223 * t246 + t232 * t247;
	t279 = t213 * t255;
	t264 = t257 * t259;
	t265 = t256 * t260;
	t230 = t254 * t250 * t256 + (t253 * t265 + t264) * t251;
	t237 = t254 * t253 - t260 * t273;
	t218 = -t230 * t246 + t237 * t247;
	t278 = t218 * t255;
	t277 = t246 * t250;
	t276 = t247 * t250;
	t275 = t247 * t255;
	t258 = cos(qJ(6));
	t274 = t247 * t258;
	t272 = t250 * t259;
	t270 = t253 * t256;
	t269 = t253 * t259;
	t266 = t256 * t257;
	t263 = t259 * t260;
	t262 = t257 * t273;
	t261 = t251 * t272;
	t236 = (-t253 * t266 + t263) * t251;
	t235 = (t253 * t264 + t265) * t251;
	t229 = t251 * t266 - t254 * t272 - t263 * t271;
	t228 = t236 * t247 + t246 * t262;
	t227 = t240 * t259 - t241 * t270;
	t226 = t240 * t256 + t241 * t269;
	t225 = t238 * t259 - t239 * t270;
	t224 = t238 * t256 + t239 * t269;
	t222 = -t240 * t269 + t241 * t256 - t249 * t261;
	t220 = -t238 * t269 + t239 * t256 + t252 * t261;
	t219 = t230 * t247 + t237 * t246;
	t217 = t218 * t258;
	t216 = t227 * t247 + t241 * t277;
	t215 = t225 * t247 + t239 * t277;
	t214 = t223 * t247 + t232 * t246;
	t212 = t221 * t247 + t231 * t246;
	t210 = t213 * t258;
	t209 = t211 * t258;
	t1 = [0, t216 * t258 + t226 * t255, -t222 * t274 + t223 * t255, t210, t210, -t214 * t255 + t222 * t258; 0, t215 * t258 + t224 * t255, -t220 * t274 + t221 * t255, t209, t209, -t212 * t255 + t220 * t258; 0, t228 * t258 + t235 * t255, -t229 * t274 + t230 * t255, t217, t217, -t219 * t255 + t229 * t258; 0, -t216 * t255 + t226 * t258, t222 * t275 + t223 * t258, -t279, -t279, -t214 * t258 - t222 * t255; 0, -t215 * t255 + t224 * t258, t220 * t275 + t221 * t258, -t280, -t280, -t212 * t258 - t220 * t255; 0, -t228 * t255 + t235 * t258, t229 * t275 + t230 * t258, -t278, -t278, -t219 * t258 - t229 * t255; 0, t227 * t246 - t241 * t276, -t222 * t246, t214, t214, 0; 0, t225 * t246 - t239 * t276, -t220 * t246, t212, t212, 0; 0, t236 * t246 - t247 * t262, -t229 * t246, t219, t219, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end