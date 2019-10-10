% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR5
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
% Datum: 2019-10-09 22:54
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR5_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR5_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:07
	% EndTime: 2019-10-09 22:54:07
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
	% StartTime: 2019-10-09 22:54:08
	% EndTime: 2019-10-09 22:54:08
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
	% StartTime: 2019-10-09 22:54:08
	% EndTime: 2019-10-09 22:54:08
	% DurationCPUTime: 0.14s
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
	% StartTime: 2019-10-09 22:54:08
	% EndTime: 2019-10-09 22:54:08
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (130->39), mult. (304->88), div. (0->0), fcn. (425->12), ass. (0->46)
	t163 = qJ(4) + pkin(13);
	t161 = sin(t163);
	t165 = sin(pkin(7));
	t189 = t161 * t165;
	t162 = cos(t163);
	t188 = t162 * t165;
	t166 = sin(pkin(6));
	t187 = t165 * t166;
	t169 = cos(pkin(6));
	t186 = t165 * t169;
	t168 = cos(pkin(7));
	t185 = t166 * t168;
	t170 = sin(qJ(3));
	t184 = t168 * t170;
	t172 = cos(qJ(3));
	t183 = t168 * t172;
	t171 = sin(qJ(2));
	t182 = t169 * t171;
	t173 = cos(qJ(2));
	t181 = t169 * t173;
	t180 = t170 * t171;
	t179 = t170 * t173;
	t178 = t171 * t172;
	t177 = t172 * t173;
	t176 = t171 * t187;
	t164 = sin(pkin(12));
	t167 = cos(pkin(12));
	t156 = -t164 * t171 + t167 * t181;
	t175 = t156 * t168 - t167 * t187;
	t158 = -t164 * t181 - t167 * t171;
	t174 = t158 * t168 + t164 * t187;
	t159 = -t164 * t182 + t167 * t173;
	t157 = t164 * t173 + t167 * t182;
	t155 = t169 * t168 - t173 * t187;
	t154 = (-t168 * t180 + t177) * t166;
	t153 = -t158 * t165 + t164 * t185;
	t152 = -t156 * t165 - t167 * t185;
	t151 = t170 * t186 + (t168 * t179 + t178) * t166;
	t150 = t172 * t186 + (t168 * t177 - t180) * t166;
	t149 = t158 * t172 - t159 * t184;
	t148 = t156 * t172 - t157 * t184;
	t147 = t159 * t172 + t174 * t170;
	t146 = -t159 * t170 + t174 * t172;
	t145 = t157 * t172 + t175 * t170;
	t144 = -t157 * t170 + t175 * t172;
	t1 = [0, t149 * t162 + t159 * t189, t146 * t162, -t147 * t161 + t153 * t162, 0, 0; 0, t148 * t162 + t157 * t189, t144 * t162, -t145 * t161 + t152 * t162, 0, 0; 0, t154 * t162 + t161 * t176, t150 * t162, -t151 * t161 + t155 * t162, 0, 0; 0, -t149 * t161 + t159 * t188, -t146 * t161, -t147 * t162 - t153 * t161, 0, 0; 0, -t148 * t161 + t157 * t188, -t144 * t161, -t145 * t162 - t152 * t161, 0, 0; 0, -t154 * t161 + t162 * t176, -t150 * t161, -t151 * t162 - t155 * t161, 0, 0; 0, t158 * t170 + t159 * t183, t147, 0, 0, 0; 0, t156 * t170 + t157 * t183, t145, 0, 0, 0; 0, (t168 * t178 + t179) * t166, t151, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:54:08
	% EndTime: 2019-10-09 22:54:09
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (288->62), mult. (691->133), div. (0->0), fcn. (952->14), ass. (0->61)
	t227 = qJ(4) + pkin(13);
	t225 = sin(t227);
	t229 = sin(pkin(7));
	t256 = t225 * t229;
	t226 = cos(t227);
	t255 = t226 * t229;
	t234 = sin(qJ(6));
	t254 = t226 * t234;
	t237 = cos(qJ(6));
	t253 = t226 * t237;
	t230 = sin(pkin(6));
	t252 = t229 * t230;
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
	t241 = t236 * t252;
	t240 = t230 * t251;
	t231 = cos(pkin(12));
	t228 = sin(pkin(12));
	t220 = -t228 * t247 + t231 * t239;
	t219 = -t228 * t246 - t231 * t236;
	t218 = t228 * t239 + t231 * t247;
	t217 = -t228 * t236 + t231 * t246;
	t216 = t233 * t232 - t239 * t252;
	t215 = (-t232 * t245 + t242) * t230;
	t214 = (t232 * t243 + t244) * t230;
	t211 = -t219 * t229 + t228 * t250;
	t210 = -t217 * t229 - t231 * t250;
	t209 = t233 * t229 * t235 + (t232 * t244 + t243) * t230;
	t208 = t230 * t245 - t233 * t251 - t242 * t250;
	t207 = t215 * t226 + t225 * t241;
	t206 = t219 * t238 - t220 * t249;
	t205 = t219 * t235 + t220 * t248;
	t204 = t217 * t238 - t218 * t249;
	t203 = t217 * t235 + t218 * t248;
	t202 = t220 * t238 + (t219 * t232 + t228 * t252) * t235;
	t201 = -t219 * t248 + t220 * t235 - t228 * t240;
	t200 = t218 * t238 + (t217 * t232 - t231 * t252) * t235;
	t199 = -t217 * t248 + t218 * t235 + t231 * t240;
	t198 = t209 * t226 + t216 * t225;
	t197 = -t209 * t225 + t216 * t226;
	t196 = t206 * t226 + t220 * t256;
	t195 = t204 * t226 + t218 * t256;
	t194 = t202 * t226 + t211 * t225;
	t193 = -t202 * t225 + t211 * t226;
	t192 = t200 * t226 + t210 * t225;
	t191 = -t200 * t225 + t210 * t226;
	t1 = [0, t196 * t237 + t205 * t234, -t201 * t253 + t202 * t234, t193 * t237, 0, -t194 * t234 + t201 * t237; 0, t195 * t237 + t203 * t234, -t199 * t253 + t200 * t234, t191 * t237, 0, -t192 * t234 + t199 * t237; 0, t207 * t237 + t214 * t234, -t208 * t253 + t209 * t234, t197 * t237, 0, -t198 * t234 + t208 * t237; 0, -t196 * t234 + t205 * t237, t201 * t254 + t202 * t237, -t193 * t234, 0, -t194 * t237 - t201 * t234; 0, -t195 * t234 + t203 * t237, t199 * t254 + t200 * t237, -t191 * t234, 0, -t192 * t237 - t199 * t234; 0, -t207 * t234 + t214 * t237, t208 * t254 + t209 * t237, -t197 * t234, 0, -t198 * t237 - t208 * t234; 0, t206 * t225 - t220 * t255, -t201 * t225, t194, 0, 0; 0, t204 * t225 - t218 * t255, -t199 * t225, t192, 0, 0; 0, t215 * t225 - t226 * t241, -t208 * t225, t198, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end