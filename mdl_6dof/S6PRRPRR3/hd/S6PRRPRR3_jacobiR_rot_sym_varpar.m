% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR3
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
% Datum: 2019-10-09 22:29
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR3_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiR_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR3_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiR_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
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
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.10s
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
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:26
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (68->24), mult. (196->58), div. (0->0), fcn. (277->12), ass. (0->27)
	t109 = sin(pkin(12));
	t111 = sin(pkin(6));
	t124 = t109 * t111;
	t113 = cos(pkin(12));
	t123 = t111 * t113;
	t115 = cos(pkin(6));
	t117 = sin(qJ(2));
	t122 = t115 * t117;
	t119 = cos(qJ(2));
	t121 = t115 * t119;
	t108 = sin(pkin(13));
	t112 = cos(pkin(13));
	t116 = sin(qJ(3));
	t118 = cos(qJ(3));
	t120 = t118 * t108 + t116 * t112;
	t104 = t116 * t108 - t118 * t112;
	t114 = cos(pkin(7));
	t110 = sin(pkin(7));
	t103 = -t109 * t122 + t113 * t119;
	t102 = -t109 * t121 - t113 * t117;
	t101 = t109 * t119 + t113 * t122;
	t100 = -t109 * t117 + t113 * t121;
	t99 = t120 * t114;
	t98 = t104 * t114;
	t97 = t120 * t110;
	t96 = t104 * t110;
	t1 = [0, -t102 * t104 - t103 * t99, -t102 * t98 - t103 * t120 - t96 * t124, 0, 0, 0; 0, -t100 * t104 - t101 * t99, -t100 * t98 - t101 * t120 + t96 * t123, 0, 0, 0; 0, (-t104 * t119 - t117 * t99) * t111, -t115 * t96 + (-t117 * t120 - t119 * t98) * t111, 0, 0, 0; 0, -t102 * t120 + t103 * t98, -t102 * t99 + t103 * t104 - t97 * t124, 0, 0, 0; 0, -t100 * t120 + t101 * t98, -t100 * t99 + t101 * t104 + t97 * t123, 0, 0, 0; 0, (t117 * t98 - t119 * t120) * t111, -t115 * t97 + (t104 * t117 - t119 * t99) * t111, 0, 0, 0; 0, t103 * t110, 0, 0, 0, 0; 0, t101 * t110, 0, 0, 0, 0; 0, t111 * t117 * t110, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:26
	% EndTime: 2019-10-09 22:29:27
	% DurationCPUTime: 0.26s
	% Computational Cost: add. (163->42), mult. (469->95), div. (0->0), fcn. (653->14), ass. (0->46)
	t196 = cos(qJ(3));
	t195 = cos(pkin(13));
	t174 = sin(pkin(12));
	t176 = sin(pkin(6));
	t194 = t174 * t176;
	t175 = sin(pkin(7));
	t193 = t175 * t176;
	t180 = sin(qJ(5));
	t192 = t175 * t180;
	t183 = cos(qJ(5));
	t191 = t175 * t183;
	t177 = cos(pkin(12));
	t190 = t176 * t177;
	t178 = cos(pkin(7));
	t189 = t176 * t178;
	t179 = cos(pkin(6));
	t182 = sin(qJ(2));
	t188 = t179 * t182;
	t184 = cos(qJ(2));
	t187 = t179 * t184;
	t186 = t182 * t193;
	t173 = sin(pkin(13));
	t181 = sin(qJ(3));
	t170 = -t196 * t173 - t181 * t195;
	t185 = -t181 * t173 + t196 * t195;
	t161 = t170 * t175;
	t163 = t170 * t178;
	t165 = -t174 * t182 + t177 * t187;
	t166 = t174 * t184 + t177 * t188;
	t150 = t161 * t190 - t165 * t163 + t166 * t185;
	t167 = -t174 * t187 - t177 * t182;
	t168 = -t174 * t188 + t177 * t184;
	t152 = -t161 * t194 - t167 * t163 + t168 * t185;
	t154 = -t179 * t161 + (-t163 * t184 + t182 * t185) * t176;
	t164 = t179 * t178 - t184 * t193;
	t162 = t185 * t178;
	t160 = t185 * t175;
	t159 = -t167 * t175 + t174 * t189;
	t158 = -t165 * t175 - t177 * t189;
	t157 = (t163 * t182 + t184 * t185) * t176;
	t156 = t168 * t163 + t167 * t185;
	t155 = t166 * t163 + t165 * t185;
	t153 = t179 * t160 + (t162 * t184 + t170 * t182) * t176;
	t151 = t160 * t194 + t167 * t162 + t168 * t170;
	t149 = -t160 * t190 + t165 * t162 + t166 * t170;
	t1 = [0, t156 * t183 + t168 * t192, t151 * t183, 0, -t152 * t180 + t159 * t183, 0; 0, t155 * t183 + t166 * t192, t149 * t183, 0, -t150 * t180 + t158 * t183, 0; 0, t157 * t183 + t180 * t186, t153 * t183, 0, -t154 * t180 + t164 * t183, 0; 0, -t156 * t180 + t168 * t191, -t151 * t180, 0, -t152 * t183 - t159 * t180, 0; 0, -t155 * t180 + t166 * t191, -t149 * t180, 0, -t150 * t183 - t158 * t180, 0; 0, -t157 * t180 + t183 * t186, -t153 * t180, 0, -t154 * t183 - t164 * t180, 0; 0, t168 * t162 - t167 * t170, t152, 0, 0, 0; 0, t166 * t162 - t165 * t170, t150, 0, 0, 0; 0, (t162 * t182 - t170 * t184) * t176, t154, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiR_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:29:27
	% EndTime: 2019-10-09 22:29:27
	% DurationCPUTime: 0.35s
	% Computational Cost: add. (369->60), mult. (1054->136), div. (0->0), fcn. (1453->16), ass. (0->62)
	t268 = cos(qJ(3));
	t267 = cos(pkin(13));
	t239 = sin(pkin(12));
	t241 = sin(pkin(6));
	t266 = t239 * t241;
	t240 = sin(pkin(7));
	t265 = t240 * t241;
	t246 = sin(qJ(5));
	t264 = t240 * t246;
	t250 = cos(qJ(5));
	t263 = t240 * t250;
	t242 = cos(pkin(12));
	t262 = t241 * t242;
	t243 = cos(pkin(7));
	t261 = t241 * t243;
	t244 = cos(pkin(6));
	t248 = sin(qJ(2));
	t260 = t244 * t248;
	t251 = cos(qJ(2));
	t259 = t244 * t251;
	t245 = sin(qJ(6));
	t258 = t245 * t250;
	t249 = cos(qJ(6));
	t257 = t249 * t250;
	t256 = t248 * t265;
	t238 = sin(pkin(13));
	t247 = sin(qJ(3));
	t235 = -t268 * t238 - t247 * t267;
	t255 = -t247 * t238 + t268 * t267;
	t226 = t235 * t240;
	t228 = t235 * t243;
	t230 = -t239 * t248 + t242 * t259;
	t231 = t239 * t251 + t242 * t260;
	t254 = t226 * t262 - t230 * t228 + t231 * t255;
	t232 = -t239 * t259 - t242 * t248;
	t233 = -t239 * t260 + t242 * t251;
	t253 = -t226 * t266 - t232 * t228 + t233 * t255;
	t252 = -t244 * t226 + (-t228 * t251 + t248 * t255) * t241;
	t229 = t244 * t243 - t251 * t265;
	t227 = t255 * t243;
	t225 = t255 * t240;
	t223 = -t232 * t240 + t239 * t261;
	t222 = -t230 * t240 - t242 * t261;
	t221 = (t228 * t248 + t251 * t255) * t241;
	t220 = (t227 * t248 - t235 * t251) * t241;
	t219 = t221 * t250 + t246 * t256;
	t218 = t233 * t228 + t232 * t255;
	t217 = t233 * t227 - t232 * t235;
	t216 = t231 * t228 + t230 * t255;
	t215 = t231 * t227 - t230 * t235;
	t213 = t244 * t225 + (t227 * t251 + t235 * t248) * t241;
	t211 = t229 * t246 + t250 * t252;
	t210 = t229 * t250 - t246 * t252;
	t208 = t225 * t266 + t232 * t227 + t233 * t235;
	t205 = -t225 * t262 + t230 * t227 + t231 * t235;
	t203 = t218 * t250 + t233 * t264;
	t202 = t216 * t250 + t231 * t264;
	t201 = t223 * t246 + t250 * t253;
	t200 = t223 * t250 - t246 * t253;
	t199 = t222 * t246 + t250 * t254;
	t198 = t222 * t250 - t246 * t254;
	t1 = [0, t203 * t249 + t217 * t245, t208 * t257 + t245 * t253, 0, t200 * t249, -t201 * t245 - t208 * t249; 0, t202 * t249 + t215 * t245, t205 * t257 + t245 * t254, 0, t198 * t249, -t199 * t245 - t205 * t249; 0, t219 * t249 + t220 * t245, t213 * t257 + t245 * t252, 0, t210 * t249, -t211 * t245 - t213 * t249; 0, -t203 * t245 + t217 * t249, -t208 * t258 + t249 * t253, 0, -t200 * t245, -t201 * t249 + t208 * t245; 0, -t202 * t245 + t215 * t249, -t205 * t258 + t249 * t254, 0, -t198 * t245, -t199 * t249 + t205 * t245; 0, -t219 * t245 + t220 * t249, -t213 * t258 + t249 * t252, 0, -t210 * t245, -t211 * t249 + t213 * t245; 0, t218 * t246 - t233 * t263, t208 * t246, 0, t201, 0; 0, t216 * t246 - t231 * t263, t205 * t246, 0, t199, 0; 0, t221 * t246 - t250 * t256, t213 * t246, 0, t211, 0;];
	JR_rot = t1;
else
	JR_rot=NaN(9,6);
end