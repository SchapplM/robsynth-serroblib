% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:04
	% EndTime: 2019-10-09 23:00:04
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
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
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
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
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
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:05
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (100->38), mult. (304->88), div. (0->0), fcn. (425->12), ass. (0->45)
	t173 = sin(pkin(7));
	t174 = sin(pkin(6));
	t199 = t173 * t174;
	t177 = cos(pkin(6));
	t198 = t173 * t177;
	t178 = sin(qJ(4));
	t197 = t173 * t178;
	t181 = cos(qJ(4));
	t196 = t173 * t181;
	t176 = cos(pkin(7));
	t195 = t174 * t176;
	t179 = sin(qJ(3));
	t194 = t176 * t179;
	t182 = cos(qJ(3));
	t193 = t176 * t182;
	t180 = sin(qJ(2));
	t192 = t177 * t180;
	t183 = cos(qJ(2));
	t191 = t177 * t183;
	t190 = t179 * t180;
	t189 = t179 * t183;
	t188 = t180 * t182;
	t187 = t182 * t183;
	t186 = t180 * t199;
	t172 = sin(pkin(12));
	t175 = cos(pkin(12));
	t167 = -t172 * t180 + t175 * t191;
	t185 = t167 * t176 - t175 * t199;
	t169 = -t172 * t191 - t175 * t180;
	t184 = t169 * t176 + t172 * t199;
	t170 = -t172 * t192 + t175 * t183;
	t168 = t172 * t183 + t175 * t192;
	t166 = t177 * t176 - t183 * t199;
	t165 = (-t176 * t190 + t187) * t174;
	t164 = -t169 * t173 + t172 * t195;
	t163 = -t167 * t173 - t175 * t195;
	t162 = t179 * t198 + (t176 * t189 + t188) * t174;
	t161 = t182 * t198 + (t176 * t187 - t190) * t174;
	t160 = t169 * t182 - t170 * t194;
	t159 = t167 * t182 - t168 * t194;
	t158 = t170 * t182 + t184 * t179;
	t157 = -t170 * t179 + t184 * t182;
	t156 = t168 * t182 + t185 * t179;
	t155 = -t168 * t179 + t185 * t182;
	t1 = [0, t169 * t179 + t170 * t193, t158, 0, 0, 0; 0, t167 * t179 + t168 * t193, t156, 0, 0, 0; 0, (t176 * t188 + t189) * t174, t162, 0, 0, 0; 0, -t160 * t181 - t170 * t197, -t157 * t181, t158 * t178 - t164 * t181, 0, 0; 0, -t159 * t181 - t168 * t197, -t155 * t181, t156 * t178 - t163 * t181, 0, 0; 0, -t165 * t181 - t178 * t186, -t161 * t181, t162 * t178 - t166 * t181, 0, 0; 0, t160 * t178 - t170 * t196, t157 * t178, t158 * t181 + t164 * t178, 0, 0; 0, t159 * t178 - t168 * t196, t155 * t178, t156 * t181 + t163 * t178, 0, 0; 0, t165 * t178 - t181 * t186, t161 * t178, t162 * t181 + t166 * t178, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:00:05
	% EndTime: 2019-10-09 23:00:06
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (228->61), mult. (691->133), div. (0->0), fcn. (952->14), ass. (0->60)
	t220 = sin(pkin(7));
	t221 = sin(pkin(6));
	t249 = t220 * t221;
	t226 = sin(qJ(4));
	t248 = t220 * t226;
	t230 = cos(qJ(4));
	t247 = t220 * t230;
	t231 = cos(qJ(3));
	t246 = t220 * t231;
	t223 = cos(pkin(7));
	t245 = t221 * t223;
	t227 = sin(qJ(3));
	t244 = t223 * t227;
	t243 = t223 * t231;
	t224 = cos(pkin(6));
	t228 = sin(qJ(2));
	t242 = t224 * t228;
	t232 = cos(qJ(2));
	t241 = t224 * t232;
	t225 = sin(qJ(6));
	t240 = t225 * t226;
	t229 = cos(qJ(6));
	t239 = t226 * t229;
	t238 = t227 * t228;
	t237 = t227 * t232;
	t236 = t228 * t231;
	t235 = t231 * t232;
	t234 = t228 * t249;
	t233 = t221 * t246;
	t222 = cos(pkin(12));
	t219 = sin(pkin(12));
	t214 = -t219 * t242 + t222 * t232;
	t213 = -t219 * t241 - t222 * t228;
	t212 = t219 * t232 + t222 * t242;
	t211 = -t219 * t228 + t222 * t241;
	t210 = t223 * t224 - t232 * t249;
	t209 = (-t223 * t238 + t235) * t221;
	t208 = (t223 * t236 + t237) * t221;
	t205 = -t213 * t220 + t219 * t245;
	t204 = -t211 * t220 - t222 * t245;
	t203 = t224 * t220 * t227 + (t223 * t237 + t236) * t221;
	t202 = t221 * t238 - t224 * t246 - t235 * t245;
	t201 = t209 * t226 - t230 * t234;
	t200 = t213 * t231 - t214 * t244;
	t199 = t213 * t227 + t214 * t243;
	t198 = t211 * t231 - t212 * t244;
	t197 = t211 * t227 + t212 * t243;
	t196 = t203 * t230 + t210 * t226;
	t195 = t203 * t226 - t210 * t230;
	t194 = t214 * t231 + (t213 * t223 + t219 * t249) * t227;
	t193 = -t213 * t243 + t214 * t227 - t219 * t233;
	t192 = t212 * t231 + (t211 * t223 - t222 * t249) * t227;
	t191 = -t211 * t243 + t212 * t227 + t222 * t233;
	t190 = t200 * t226 - t214 * t247;
	t189 = t198 * t226 - t212 * t247;
	t188 = t194 * t230 + t205 * t226;
	t187 = t194 * t226 - t205 * t230;
	t186 = t192 * t230 + t204 * t226;
	t185 = t192 * t226 - t204 * t230;
	t1 = [0, t190 * t225 + t199 * t229, -t193 * t240 + t194 * t229, t188 * t225, 0, t187 * t229 - t193 * t225; 0, t189 * t225 + t197 * t229, -t191 * t240 + t192 * t229, t186 * t225, 0, t185 * t229 - t191 * t225; 0, t201 * t225 + t208 * t229, -t202 * t240 + t203 * t229, t196 * t225, 0, t195 * t229 - t202 * t225; 0, t190 * t229 - t199 * t225, -t193 * t239 - t194 * t225, t188 * t229, 0, -t187 * t225 - t193 * t229; 0, t189 * t229 - t197 * t225, -t191 * t239 - t192 * t225, t186 * t229, 0, -t185 * t225 - t191 * t229; 0, t201 * t229 - t208 * t225, -t202 * t239 - t203 * t225, t196 * t229, 0, -t195 * t225 - t202 * t229; 0, t200 * t230 + t214 * t248, -t193 * t230, -t187, 0, 0; 0, t198 * t230 + t212 * t248, -t191 * t230, -t185, 0, 0; 0, t209 * t230 + t226 * t234, -t202 * t230, -t195, 0, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end