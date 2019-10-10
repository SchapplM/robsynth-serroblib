% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:35
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR6_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR6_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
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
	% StartTime: 2019-10-09 22:35:19
	% EndTime: 2019-10-09 22:35:19
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
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:20
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (64->29), mult. (200->69), div. (0->0), fcn. (279->12), ass. (0->38)
	t139 = sin(pkin(13));
	t141 = sin(pkin(7));
	t165 = t139 * t141;
	t142 = sin(pkin(6));
	t164 = t141 * t142;
	t143 = cos(pkin(13));
	t163 = t141 * t143;
	t146 = cos(pkin(6));
	t162 = t141 * t146;
	t145 = cos(pkin(7));
	t147 = sin(qJ(3));
	t161 = t145 * t147;
	t149 = cos(qJ(3));
	t160 = t145 * t149;
	t148 = sin(qJ(2));
	t159 = t146 * t148;
	t150 = cos(qJ(2));
	t158 = t146 * t150;
	t157 = t147 * t148;
	t156 = t147 * t150;
	t155 = t148 * t149;
	t154 = t149 * t150;
	t153 = t148 * t164;
	t140 = sin(pkin(12));
	t144 = cos(pkin(12));
	t134 = -t140 * t148 + t144 * t158;
	t152 = t134 * t145 - t144 * t164;
	t136 = -t140 * t158 - t144 * t148;
	t151 = t136 * t145 + t140 * t164;
	t137 = -t140 * t159 + t144 * t150;
	t135 = t140 * t150 + t144 * t159;
	t133 = (-t145 * t157 + t154) * t142;
	t132 = t149 * t162 + (t145 * t154 - t157) * t142;
	t131 = t136 * t149 - t137 * t161;
	t130 = t134 * t149 - t135 * t161;
	t129 = -t137 * t147 + t151 * t149;
	t128 = -t135 * t147 + t152 * t149;
	t1 = [0, t131 * t143 + t137 * t165, t129 * t143, 0, 0, 0; 0, t130 * t143 + t135 * t165, t128 * t143, 0, 0, 0; 0, t133 * t143 + t139 * t153, t132 * t143, 0, 0, 0; 0, -t131 * t139 + t137 * t163, -t129 * t139, 0, 0, 0; 0, -t130 * t139 + t135 * t163, -t128 * t139, 0, 0, 0; 0, -t133 * t139 + t143 * t153, -t132 * t139, 0, 0, 0; 0, t136 * t147 + t137 * t160, t137 * t149 + t151 * t147, 0, 0, 0; 0, t134 * t147 + t135 * t160, t135 * t149 + t152 * t147, 0, 0, 0; 0, (t145 * t155 + t156) * t142, t147 * t162 + (t145 * t156 + t155) * t142, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:20
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (130->39), mult. (304->88), div. (0->0), fcn. (425->12), ass. (0->46)
	t160 = pkin(13) + qJ(5);
	t158 = sin(t160);
	t162 = sin(pkin(7));
	t186 = t158 * t162;
	t159 = cos(t160);
	t185 = t159 * t162;
	t163 = sin(pkin(6));
	t184 = t162 * t163;
	t166 = cos(pkin(6));
	t183 = t162 * t166;
	t165 = cos(pkin(7));
	t182 = t163 * t165;
	t167 = sin(qJ(3));
	t181 = t165 * t167;
	t169 = cos(qJ(3));
	t180 = t165 * t169;
	t168 = sin(qJ(2));
	t179 = t166 * t168;
	t170 = cos(qJ(2));
	t178 = t166 * t170;
	t177 = t167 * t168;
	t176 = t167 * t170;
	t175 = t168 * t169;
	t174 = t169 * t170;
	t173 = t168 * t184;
	t161 = sin(pkin(12));
	t164 = cos(pkin(12));
	t153 = -t161 * t168 + t164 * t178;
	t172 = t153 * t165 - t164 * t184;
	t155 = -t161 * t178 - t164 * t168;
	t171 = t155 * t165 + t161 * t184;
	t156 = -t161 * t179 + t164 * t170;
	t154 = t161 * t170 + t164 * t179;
	t152 = t166 * t165 - t170 * t184;
	t151 = (-t165 * t177 + t174) * t163;
	t150 = -t155 * t162 + t161 * t182;
	t149 = -t153 * t162 - t164 * t182;
	t148 = t167 * t183 + (t165 * t176 + t175) * t163;
	t147 = t169 * t183 + (t165 * t174 - t177) * t163;
	t146 = t155 * t169 - t156 * t181;
	t145 = t153 * t169 - t154 * t181;
	t144 = t156 * t169 + t171 * t167;
	t143 = -t156 * t167 + t171 * t169;
	t142 = t154 * t169 + t172 * t167;
	t141 = -t154 * t167 + t172 * t169;
	t1 = [0, t146 * t159 + t156 * t186, t143 * t159, 0, -t144 * t158 + t150 * t159, 0; 0, t145 * t159 + t154 * t186, t141 * t159, 0, -t142 * t158 + t149 * t159, 0; 0, t151 * t159 + t158 * t173, t147 * t159, 0, -t148 * t158 + t152 * t159, 0; 0, -t146 * t158 + t156 * t185, -t143 * t158, 0, -t144 * t159 - t150 * t158, 0; 0, -t145 * t158 + t154 * t185, -t141 * t158, 0, -t142 * t159 - t149 * t158, 0; 0, -t151 * t158 + t159 * t173, -t147 * t158, 0, -t148 * t159 - t152 * t158, 0; 0, t155 * t167 + t156 * t180, t144, 0, 0, 0; 0, t153 * t167 + t154 * t180, t142, 0, 0, 0; 0, (t165 * t175 + t176) * t163, t148, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:35:20
	% EndTime: 2019-10-09 22:35:20
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (288->62), mult. (691->133), div. (0->0), fcn. (952->14), ass. (0->61)
	t226 = pkin(13) + qJ(5);
	t224 = sin(t226);
	t228 = sin(pkin(7));
	t255 = t224 * t228;
	t225 = cos(t226);
	t254 = t225 * t228;
	t233 = sin(qJ(6));
	t253 = t225 * t233;
	t236 = cos(qJ(6));
	t252 = t225 * t236;
	t229 = sin(pkin(6));
	t251 = t228 * t229;
	t237 = cos(qJ(3));
	t250 = t228 * t237;
	t231 = cos(pkin(7));
	t249 = t229 * t231;
	t234 = sin(qJ(3));
	t248 = t231 * t234;
	t247 = t231 * t237;
	t232 = cos(pkin(6));
	t235 = sin(qJ(2));
	t246 = t232 * t235;
	t238 = cos(qJ(2));
	t245 = t232 * t238;
	t244 = t234 * t235;
	t243 = t234 * t238;
	t242 = t235 * t237;
	t241 = t237 * t238;
	t240 = t235 * t251;
	t239 = t229 * t250;
	t230 = cos(pkin(12));
	t227 = sin(pkin(12));
	t219 = -t227 * t246 + t230 * t238;
	t218 = -t227 * t245 - t230 * t235;
	t217 = t227 * t238 + t230 * t246;
	t216 = -t227 * t235 + t230 * t245;
	t215 = t232 * t231 - t238 * t251;
	t214 = (-t231 * t244 + t241) * t229;
	t213 = (t231 * t242 + t243) * t229;
	t210 = -t218 * t228 + t227 * t249;
	t209 = -t216 * t228 - t230 * t249;
	t208 = t232 * t228 * t234 + (t231 * t243 + t242) * t229;
	t207 = t229 * t244 - t232 * t250 - t241 * t249;
	t206 = t214 * t225 + t224 * t240;
	t205 = t218 * t237 - t219 * t248;
	t204 = t218 * t234 + t219 * t247;
	t203 = t216 * t237 - t217 * t248;
	t202 = t216 * t234 + t217 * t247;
	t201 = t219 * t237 + (t218 * t231 + t227 * t251) * t234;
	t200 = -t218 * t247 + t219 * t234 - t227 * t239;
	t199 = t217 * t237 + (t216 * t231 - t230 * t251) * t234;
	t198 = -t216 * t247 + t217 * t234 + t230 * t239;
	t197 = t208 * t225 + t215 * t224;
	t196 = -t208 * t224 + t215 * t225;
	t195 = t205 * t225 + t219 * t255;
	t194 = t203 * t225 + t217 * t255;
	t193 = t201 * t225 + t210 * t224;
	t192 = -t201 * t224 + t210 * t225;
	t191 = t199 * t225 + t209 * t224;
	t190 = -t199 * t224 + t209 * t225;
	t1 = [0, t195 * t236 + t204 * t233, -t200 * t252 + t201 * t233, 0, t192 * t236, -t193 * t233 + t200 * t236; 0, t194 * t236 + t202 * t233, -t198 * t252 + t199 * t233, 0, t190 * t236, -t191 * t233 + t198 * t236; 0, t206 * t236 + t213 * t233, -t207 * t252 + t208 * t233, 0, t196 * t236, -t197 * t233 + t207 * t236; 0, -t195 * t233 + t204 * t236, t200 * t253 + t201 * t236, 0, -t192 * t233, -t193 * t236 - t200 * t233; 0, -t194 * t233 + t202 * t236, t198 * t253 + t199 * t236, 0, -t190 * t233, -t191 * t236 - t198 * t233; 0, -t206 * t233 + t213 * t236, t207 * t253 + t208 * t236, 0, -t196 * t233, -t197 * t236 - t207 * t233; 0, t205 * t224 - t219 * t254, -t200 * t224, 0, t193, 0; 0, t203 * t224 - t217 * t254, -t198 * t224, 0, t191, 0; 0, t214 * t224 - t225 * t240, -t207 * t224, 0, t197, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end