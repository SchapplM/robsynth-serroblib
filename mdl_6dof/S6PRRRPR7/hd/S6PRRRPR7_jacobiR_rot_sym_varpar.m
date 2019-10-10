% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR7_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR7_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(6));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(12));
	t20 = sin(pkin(6));
	t19 = sin(pkin(12));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0, 0; 0, t20 * t24, 0, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0, 0; 0, -t20 * t23, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
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
	t80 = sin(pkin(12));
	t83 = cos(pkin(12));
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
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:04
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
	t151 = sin(pkin(12));
	t154 = cos(pkin(12));
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
	% StartTime: 2019-10-09 22:58:04
	% EndTime: 2019-10-09 22:58:05
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (169->50), mult. (513->117), div. (0->0), fcn. (706->14), ass. (0->58)
	t195 = sin(pkin(13));
	t206 = cos(qJ(4));
	t226 = t195 * t206;
	t197 = sin(pkin(7));
	t198 = sin(pkin(6));
	t225 = t197 * t198;
	t202 = cos(pkin(6));
	t224 = t197 * t202;
	t203 = sin(qJ(4));
	t223 = t197 * t203;
	t222 = t197 * t206;
	t201 = cos(pkin(7));
	t221 = t198 * t201;
	t199 = cos(pkin(13));
	t220 = t199 * t206;
	t204 = sin(qJ(3));
	t219 = t201 * t204;
	t207 = cos(qJ(3));
	t218 = t201 * t207;
	t205 = sin(qJ(2));
	t217 = t202 * t205;
	t208 = cos(qJ(2));
	t216 = t202 * t208;
	t215 = t204 * t205;
	t214 = t204 * t208;
	t213 = t205 * t207;
	t212 = t207 * t208;
	t211 = t205 * t225;
	t196 = sin(pkin(12));
	t200 = cos(pkin(12));
	t190 = -t196 * t205 + t200 * t216;
	t210 = t190 * t201 - t200 * t225;
	t192 = -t196 * t216 - t200 * t205;
	t209 = t192 * t201 + t196 * t225;
	t193 = -t196 * t217 + t200 * t208;
	t191 = t196 * t208 + t200 * t217;
	t189 = t202 * t201 - t208 * t225;
	t188 = (-t201 * t215 + t212) * t198;
	t187 = (t201 * t213 + t214) * t198;
	t186 = -t192 * t197 + t196 * t221;
	t185 = -t190 * t197 - t200 * t221;
	t184 = t204 * t224 + (t201 * t214 + t213) * t198;
	t183 = t207 * t224 + (t201 * t212 - t215) * t198;
	t182 = t188 * t206 + t203 * t211;
	t181 = t192 * t207 - t193 * t219;
	t180 = t192 * t204 + t193 * t218;
	t179 = t190 * t207 - t191 * t219;
	t178 = t190 * t204 + t191 * t218;
	t177 = -t184 * t203 + t189 * t206;
	t176 = t193 * t207 + t209 * t204;
	t175 = -t193 * t204 + t209 * t207;
	t174 = t191 * t207 + t210 * t204;
	t173 = -t191 * t204 + t210 * t207;
	t172 = t181 * t206 + t193 * t223;
	t171 = t179 * t206 + t191 * t223;
	t170 = -t176 * t203 + t186 * t206;
	t169 = -t174 * t203 + t185 * t206;
	t1 = [0, t172 * t199 + t180 * t195, t175 * t220 + t176 * t195, t170 * t199, 0, 0; 0, t171 * t199 + t178 * t195, t173 * t220 + t174 * t195, t169 * t199, 0, 0; 0, t182 * t199 + t187 * t195, t183 * t220 + t184 * t195, t177 * t199, 0, 0; 0, -t172 * t195 + t180 * t199, -t175 * t226 + t176 * t199, -t170 * t195, 0, 0; 0, -t171 * t195 + t178 * t199, -t173 * t226 + t174 * t199, -t169 * t195, 0, 0; 0, -t182 * t195 + t187 * t199, -t183 * t226 + t184 * t199, -t177 * t195, 0, 0; 0, t181 * t203 - t193 * t222, t175 * t203, t176 * t206 + t186 * t203, 0, 0; 0, t179 * t203 - t191 * t222, t173 * t203, t174 * t206 + t185 * t203, 0, 0; 0, t188 * t203 - t206 * t211, t183 * t203, t184 * t206 + t189 * t203, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:58:05
	% EndTime: 2019-10-09 22:58:05
	% DurationCPUTime: 0.23s
	% Computational Cost: add. (273->62), mult. (691->133), div. (0->0), fcn. (952->14), ass. (0->61)
	t227 = pkin(13) + qJ(6);
	t225 = sin(t227);
	t237 = cos(qJ(4));
	t256 = t225 * t237;
	t226 = cos(t227);
	t255 = t226 * t237;
	t229 = sin(pkin(7));
	t230 = sin(pkin(6));
	t254 = t229 * t230;
	t234 = sin(qJ(4));
	t253 = t229 * t234;
	t252 = t229 * t237;
	t238 = cos(qJ(3));
	t251 = t229 * t238;
	t232 = cos(pkin(7));
	t250 = t230 * t232;
	t235 = sin(qJ(3));
	t249 = t232 * t235;
	t248 = t232 * t238;
	t233 = cos(pkin(6));
	t236 = sin(qJ(2));
	t247 = t233 * t236;
	t239 = cos(qJ(2));
	t246 = t233 * t239;
	t245 = t235 * t236;
	t244 = t235 * t239;
	t243 = t236 * t238;
	t242 = t238 * t239;
	t241 = t236 * t254;
	t240 = t230 * t251;
	t231 = cos(pkin(12));
	t228 = sin(pkin(12));
	t220 = -t228 * t247 + t231 * t239;
	t219 = -t228 * t246 - t231 * t236;
	t218 = t228 * t239 + t231 * t247;
	t217 = -t228 * t236 + t231 * t246;
	t216 = t233 * t232 - t239 * t254;
	t215 = (-t232 * t245 + t242) * t230;
	t214 = (t232 * t243 + t244) * t230;
	t211 = -t219 * t229 + t228 * t250;
	t210 = -t217 * t229 - t231 * t250;
	t209 = t233 * t229 * t235 + (t232 * t244 + t243) * t230;
	t208 = t230 * t245 - t233 * t251 - t242 * t250;
	t207 = t215 * t237 + t234 * t241;
	t206 = t219 * t238 - t220 * t249;
	t205 = t219 * t235 + t220 * t248;
	t204 = t217 * t238 - t218 * t249;
	t203 = t217 * t235 + t218 * t248;
	t202 = t209 * t237 + t216 * t234;
	t201 = -t209 * t234 + t216 * t237;
	t200 = t220 * t238 + (t219 * t232 + t228 * t254) * t235;
	t199 = -t219 * t248 + t220 * t235 - t228 * t240;
	t198 = t218 * t238 + (t217 * t232 - t231 * t254) * t235;
	t197 = -t217 * t248 + t218 * t235 + t231 * t240;
	t196 = t206 * t237 + t220 * t253;
	t195 = t204 * t237 + t218 * t253;
	t194 = t200 * t237 + t211 * t234;
	t193 = -t200 * t234 + t211 * t237;
	t192 = t198 * t237 + t210 * t234;
	t191 = -t198 * t234 + t210 * t237;
	t1 = [0, t196 * t226 + t205 * t225, -t199 * t255 + t200 * t225, t193 * t226, 0, -t194 * t225 + t199 * t226; 0, t195 * t226 + t203 * t225, -t197 * t255 + t198 * t225, t191 * t226, 0, -t192 * t225 + t197 * t226; 0, t207 * t226 + t214 * t225, -t208 * t255 + t209 * t225, t201 * t226, 0, -t202 * t225 + t208 * t226; 0, -t196 * t225 + t205 * t226, t199 * t256 + t200 * t226, -t193 * t225, 0, -t194 * t226 - t199 * t225; 0, -t195 * t225 + t203 * t226, t197 * t256 + t198 * t226, -t191 * t225, 0, -t192 * t226 - t197 * t225; 0, -t207 * t225 + t214 * t226, t208 * t256 + t209 * t226, -t201 * t225, 0, -t202 * t226 - t208 * t225; 0, t206 * t234 - t220 * t252, -t199 * t234, t194, 0, 0; 0, t204 * t234 - t218 * t252, -t197 * t234, t192, 0, 0; 0, t215 * t234 - t237 * t241, -t208 * t234, t202, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end