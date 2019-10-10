% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR8_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR8_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiR_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:09
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
	% StartTime: 2019-10-09 22:39:09
	% EndTime: 2019-10-09 22:39:10
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
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (38->20), mult. (118->51), div. (0->0), fcn. (169->10), ass. (0->27)
	t116 = sin(pkin(7));
	t117 = sin(pkin(6));
	t136 = t116 * t117;
	t120 = cos(pkin(6));
	t135 = t116 * t120;
	t119 = cos(pkin(7));
	t121 = sin(qJ(3));
	t134 = t119 * t121;
	t123 = cos(qJ(3));
	t133 = t119 * t123;
	t122 = sin(qJ(2));
	t132 = t120 * t122;
	t124 = cos(qJ(2));
	t131 = t120 * t124;
	t130 = t121 * t122;
	t129 = t121 * t124;
	t128 = t122 * t123;
	t127 = t123 * t124;
	t115 = sin(pkin(12));
	t118 = cos(pkin(12));
	t110 = -t115 * t122 + t118 * t131;
	t126 = -t110 * t119 + t118 * t136;
	t112 = -t115 * t131 - t118 * t122;
	t125 = t112 * t119 + t115 * t136;
	t113 = -t115 * t132 + t118 * t124;
	t111 = t115 * t124 + t118 * t132;
	t1 = [0, t113 * t116, 0, 0, 0, 0; 0, t111 * t116, 0, 0, 0, 0; 0, t122 * t136, 0, 0, 0, 0; 0, -t112 * t123 + t113 * t134, t113 * t121 - t125 * t123, 0, 0, 0; 0, -t110 * t123 + t111 * t134, t111 * t121 + t126 * t123, 0, 0, 0; 0, (t119 * t130 - t127) * t117, -t123 * t135 + (-t119 * t127 + t130) * t117, 0, 0, 0; 0, t112 * t121 + t113 * t133, t113 * t123 + t125 * t121, 0, 0, 0; 0, t110 * t121 + t111 * t133, t111 * t123 - t126 * t121, 0, 0, 0; 0, (t119 * t128 + t129) * t117, t121 * t135 + (t119 * t129 + t128) * t117, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.19s
	% Computational Cost: add. (97->38), mult. (304->88), div. (0->0), fcn. (425->12), ass. (0->45)
	t153 = sin(pkin(7));
	t154 = sin(pkin(6));
	t179 = t153 * t154;
	t157 = cos(pkin(6));
	t178 = t153 * t157;
	t158 = sin(qJ(5));
	t177 = t153 * t158;
	t161 = cos(qJ(5));
	t176 = t153 * t161;
	t156 = cos(pkin(7));
	t175 = t154 * t156;
	t159 = sin(qJ(3));
	t174 = t156 * t159;
	t162 = cos(qJ(3));
	t173 = t156 * t162;
	t160 = sin(qJ(2));
	t172 = t157 * t160;
	t163 = cos(qJ(2));
	t171 = t157 * t163;
	t170 = t159 * t160;
	t169 = t159 * t163;
	t168 = t160 * t162;
	t167 = t162 * t163;
	t166 = t160 * t179;
	t152 = sin(pkin(12));
	t155 = cos(pkin(12));
	t147 = -t152 * t160 + t155 * t171;
	t165 = -t147 * t156 + t155 * t179;
	t149 = -t152 * t171 - t155 * t160;
	t164 = t149 * t156 + t152 * t179;
	t150 = -t152 * t172 + t155 * t163;
	t148 = t152 * t163 + t155 * t172;
	t146 = t156 * t157 - t163 * t179;
	t145 = (t156 * t168 + t169) * t154;
	t144 = -t149 * t153 + t152 * t175;
	t143 = -t147 * t153 - t155 * t175;
	t142 = t159 * t178 + (t156 * t169 + t168) * t154;
	t141 = -t162 * t178 + (-t156 * t167 + t170) * t154;
	t140 = t149 * t159 + t150 * t173;
	t139 = t147 * t159 + t148 * t173;
	t138 = t150 * t162 + t159 * t164;
	t137 = t150 * t159 - t162 * t164;
	t136 = t148 * t162 - t159 * t165;
	t135 = t148 * t159 + t162 * t165;
	t1 = [0, t140 * t158 + t150 * t176, t138 * t158, 0, t137 * t161 - t144 * t158, 0; 0, t139 * t158 + t148 * t176, t136 * t158, 0, t135 * t161 - t143 * t158, 0; 0, t145 * t158 + t161 * t166, t142 * t158, 0, t141 * t161 - t146 * t158, 0; 0, t140 * t161 - t150 * t177, t138 * t161, 0, -t137 * t158 - t144 * t161, 0; 0, t139 * t161 - t148 * t177, t136 * t161, 0, -t135 * t158 - t143 * t161, 0; 0, t145 * t161 - t158 * t166, t142 * t161, 0, -t141 * t158 - t146 * t161, 0; 0, t149 * t162 - t150 * t174, -t137, 0, 0, 0; 0, t147 * t162 - t148 * t174, -t135, 0, 0, 0; 0, (-t156 * t170 + t167) * t154, -t141, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:39:10
	% EndTime: 2019-10-09 22:39:10
	% DurationCPUTime: 0.24s
	% Computational Cost: add. (234->61), mult. (691->133), div. (0->0), fcn. (952->14), ass. (0->60)
	t218 = sin(pkin(7));
	t219 = sin(pkin(6));
	t247 = t218 * t219;
	t224 = sin(qJ(5));
	t246 = t218 * t224;
	t228 = cos(qJ(5));
	t245 = t218 * t228;
	t229 = cos(qJ(3));
	t244 = t218 * t229;
	t221 = cos(pkin(7));
	t243 = t219 * t221;
	t225 = sin(qJ(3));
	t242 = t221 * t225;
	t241 = t221 * t229;
	t222 = cos(pkin(6));
	t226 = sin(qJ(2));
	t240 = t222 * t226;
	t230 = cos(qJ(2));
	t239 = t222 * t230;
	t223 = sin(qJ(6));
	t238 = t223 * t224;
	t227 = cos(qJ(6));
	t237 = t224 * t227;
	t236 = t225 * t226;
	t235 = t225 * t230;
	t234 = t226 * t229;
	t233 = t229 * t230;
	t232 = t226 * t247;
	t231 = t219 * t244;
	t220 = cos(pkin(12));
	t217 = sin(pkin(12));
	t212 = -t217 * t240 + t220 * t230;
	t211 = -t217 * t239 - t220 * t226;
	t210 = t217 * t230 + t220 * t240;
	t209 = -t217 * t226 + t220 * t239;
	t208 = t222 * t221 - t230 * t247;
	t207 = (-t221 * t236 + t233) * t219;
	t206 = (t221 * t234 + t235) * t219;
	t203 = -t211 * t218 + t217 * t243;
	t202 = -t209 * t218 - t220 * t243;
	t201 = t222 * t218 * t225 + (t221 * t235 + t234) * t219;
	t200 = t219 * t236 - t222 * t244 - t233 * t243;
	t199 = t206 * t224 + t228 * t232;
	t198 = t211 * t229 - t212 * t242;
	t197 = t211 * t225 + t212 * t241;
	t196 = t209 * t229 - t210 * t242;
	t195 = t209 * t225 + t210 * t241;
	t194 = t200 * t224 + t208 * t228;
	t193 = t200 * t228 - t208 * t224;
	t192 = t212 * t229 + (t211 * t221 + t217 * t247) * t225;
	t191 = -t211 * t241 + t212 * t225 - t217 * t231;
	t190 = t210 * t229 + (t209 * t221 - t220 * t247) * t225;
	t189 = -t209 * t241 + t210 * t225 + t220 * t231;
	t188 = t197 * t224 + t212 * t245;
	t187 = t195 * t224 + t210 * t245;
	t186 = t191 * t224 + t203 * t228;
	t185 = t191 * t228 - t203 * t224;
	t184 = t189 * t224 + t202 * t228;
	t183 = t189 * t228 - t202 * t224;
	t1 = [0, t188 * t227 + t198 * t223, -t191 * t223 + t192 * t237, 0, t185 * t227, -t186 * t223 + t192 * t227; 0, t187 * t227 + t196 * t223, -t189 * t223 + t190 * t237, 0, t183 * t227, -t184 * t223 + t190 * t227; 0, t199 * t227 + t207 * t223, -t200 * t223 + t201 * t237, 0, t193 * t227, -t194 * t223 + t201 * t227; 0, -t188 * t223 + t198 * t227, -t191 * t227 - t192 * t238, 0, -t185 * t223, -t186 * t227 - t192 * t223; 0, -t187 * t223 + t196 * t227, -t189 * t227 - t190 * t238, 0, -t183 * t223, -t184 * t227 - t190 * t223; 0, -t199 * t223 + t207 * t227, -t200 * t227 - t201 * t238, 0, -t193 * t223, -t194 * t227 - t201 * t223; 0, -t197 * t228 + t212 * t246, -t192 * t228, 0, t186, 0; 0, -t195 * t228 + t210 * t246, -t190 * t228, 0, t184, 0; 0, -t206 * t228 + t224 * t232, -t201 * t228, 0, t194, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end