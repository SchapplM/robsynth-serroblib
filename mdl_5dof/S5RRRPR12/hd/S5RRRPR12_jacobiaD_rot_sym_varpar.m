% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRRPR12_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRRPR12_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR12_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:42:20
	% EndTime: 2019-12-31 21:42:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:42:20
	% EndTime: 2019-12-31 21:42:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:42:20
	% EndTime: 2019-12-31 21:42:21
	% DurationCPUTime: 0.29s
	% Computational Cost: add. (215->39), mult. (853->106), div. (126->12), fcn. (1047->9), ass. (0->54)
	t99 = sin(pkin(5));
	t93 = t99 ^ 2;
	t100 = cos(pkin(5));
	t95 = 0.1e1 / t100 ^ 2;
	t104 = cos(qJ(1));
	t98 = t104 ^ 2;
	t89 = t98 * t93 * t95 + 0.1e1;
	t102 = sin(qJ(1));
	t97 = t102 ^ 2;
	t126 = 0.1e1 / t89 ^ 2 * t97;
	t131 = t126 * t95;
	t122 = t104 * t99;
	t88 = atan2(t122, t100);
	t84 = sin(t88);
	t85 = cos(t88);
	t72 = t85 * t100 + t84 * t122;
	t67 = 0.1e1 / t72;
	t103 = cos(qJ(2));
	t118 = t104 * t103;
	t101 = sin(qJ(2));
	t121 = t102 * t101;
	t113 = t100 * t121 - t118;
	t77 = 0.1e1 / t113;
	t94 = 0.1e1 / t100;
	t68 = 0.1e1 / t72 ^ 2;
	t78 = 0.1e1 / t113 ^ 2;
	t119 = t104 * t101;
	t120 = t102 * t103;
	t81 = -t100 * t119 - t120;
	t82 = t100 * t120 + t119;
	t71 = t81 * qJD(1) - t82 * qJD(2);
	t128 = t71 * t77 * t78;
	t115 = t100 * t118;
	t70 = -qJD(1) * t115 - qJD(2) * t118 + (qJD(2) * t100 + qJD(1)) * t121;
	t129 = t70 * t78;
	t76 = t82 ^ 2;
	t75 = t76 * t78 + 0.1e1;
	t130 = (t76 * t128 - t82 * t129) / t75 ^ 2;
	t127 = t81 * t82;
	t125 = t93 * t94;
	t124 = t102 * t68;
	t123 = t104 * t68;
	t117 = qJD(1) * t104;
	t86 = 0.1e1 / t89;
	t116 = (t86 - 0.1e1) * t99;
	t114 = -0.2e1 * t94 * t131;
	t80 = t115 - t121;
	t63 = (-t104 * t85 * t86 * t125 + t84 * t116) * t102;
	t92 = t99 * t93;
	t73 = 0.1e1 / t75;
	t69 = t67 * t68;
	t66 = t97 * t93 * t68 + 0.1e1;
	t62 = qJD(1) * t63;
	t1 = [(-t86 * t94 * t99 + t92 * t114) * t117, 0, 0, 0, 0; (0.2e1 * (-t104 * t67 + t63 * t124) / t66 ^ 2 * (-t62 * t69 * t97 + t117 * t124) * t93 + ((0.2e1 * t102 * t63 * t69 - t123) * t62 + (-t63 * t123 + (-t67 + (-t92 * t131 - t116) * t84 * t123 - (t93 ^ 2 * t98 * t114 + (-t126 + (0.2e1 * t97 - t98) * t86) * t125) * t68 * t85) * t102) * qJD(1)) / t66) * t99, 0, 0, 0, 0; 0.2e1 * (t78 * t127 + t77 * t80) * t130 + (-(-t82 * qJD(1) + t81 * qJD(2)) * t77 - 0.2e1 * t127 * t128 + (-t80 * t71 - (t113 * qJD(1) - t80 * qJD(2)) * t82 + t81 * t70) * t78) * t73, -0.2e1 * t130 + 0.2e1 * (-t73 * t129 + (t73 * t128 - t78 * t130) * t82) * t82, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:42:21
	% EndTime: 2019-12-31 21:42:21
	% DurationCPUTime: 0.76s
	% Computational Cost: add. (1479->91), mult. (4303->201), div. (668->14), fcn. (5516->11), ass. (0->91)
	t171 = sin(qJ(2));
	t172 = sin(qJ(1));
	t225 = cos(pkin(5));
	t195 = t172 * t225;
	t193 = t171 * t195;
	t174 = cos(qJ(2));
	t175 = cos(qJ(1));
	t209 = t175 * t174;
	t157 = -t193 + t209;
	t170 = sin(qJ(3));
	t173 = cos(qJ(3));
	t169 = sin(pkin(5));
	t213 = t169 * t172;
	t185 = -t157 * t170 + t173 * t213;
	t227 = t185 * qJD(3);
	t194 = t175 * t225;
	t192 = t174 * t194;
	t210 = t172 * t171;
	t153 = -t192 + t210;
	t212 = t169 * t174;
	t147 = atan2(-t153, -t212);
	t145 = sin(t147);
	t146 = cos(t147);
	t151 = t153 ^ 2;
	t165 = 0.1e1 / t169 ^ 2;
	t167 = 0.1e1 / t174 ^ 2;
	t150 = t151 * t165 * t167 + 0.1e1;
	t148 = 0.1e1 / t150;
	t164 = 0.1e1 / t169;
	t166 = 0.1e1 / t174;
	t199 = t153 * t164 * t166;
	t226 = (t146 * t199 - t145) * t148 + t145;
	t129 = -t145 * t153 - t146 * t212;
	t126 = 0.1e1 / t129;
	t144 = t157 * t173 + t170 * t213;
	t138 = 0.1e1 / t144;
	t127 = 0.1e1 / t129 ^ 2;
	t139 = 0.1e1 / t144 ^ 2;
	t182 = -t171 * t194 - t172 * t174;
	t183 = -t175 * t171 - t174 * t195;
	t135 = -t183 * qJD(1) - t182 * qJD(2);
	t207 = qJD(2) * t171;
	t196 = t167 * t207;
	t184 = t135 * t166 + t153 * t196;
	t215 = t148 * t164;
	t118 = t184 * t215;
	t188 = t145 * t212 - t146 * t153;
	t200 = t146 * t169 * t171;
	t114 = qJD(2) * t200 + t188 * t118 - t145 * t135;
	t224 = t114 * t126 * t127;
	t214 = t167 * t171;
	t187 = t153 * t214 - t166 * t182;
	t119 = t187 * t215;
	t115 = t188 * t119 + t145 * t182 + t200;
	t223 = t115 * t183;
	t134 = t182 * qJD(1) + t183 * qJD(2);
	t208 = qJD(1) * t169;
	t197 = t175 * t208;
	t124 = t144 * qJD(3) + t134 * t170 - t173 * t197;
	t137 = t185 ^ 2;
	t132 = t137 * t139 + 0.1e1;
	t218 = t139 * t185;
	t125 = t134 * t173 + t170 * t197 + t227;
	t220 = t125 * t138 * t139;
	t222 = (-t124 * t218 - t137 * t220) / t132 ^ 2;
	t168 = t166 * t167;
	t221 = (t135 * t153 * t167 + t151 * t168 * t207) * t165 / t150 ^ 2;
	t191 = qJD(2) * t225 + qJD(1);
	t206 = qJD(2) * t174;
	t133 = -qJD(1) * t192 - t175 * t206 + t191 * t210;
	t219 = t133 * t127;
	t217 = t145 * t183;
	t216 = t146 * t183;
	t211 = t169 * t175;
	t152 = t183 ^ 2;
	t122 = t127 * t152 + 0.1e1;
	t205 = 0.2e1 * (-t152 * t224 + t183 * t219) / t122 ^ 2;
	t204 = 0.2e1 * t224;
	t203 = 0.2e1 * t222;
	t202 = -0.2e1 * t221;
	t201 = t185 * t220;
	t198 = t172 * t208;
	t189 = t170 * t138 + t173 * t218;
	t186 = -t170 * t182 + t173 * t211;
	t142 = t170 * t211 + t173 * t182;
	t136 = -qJD(1) * t193 - t172 * t207 + t191 * t209;
	t130 = 0.1e1 / t132;
	t120 = 0.1e1 / t122;
	t117 = t226 * t183;
	t113 = (t187 * t202 + (t135 * t214 + t136 * t166 + (-t182 * t214 + (0.2e1 * t168 * t171 ^ 2 + t166) * t153) * qJD(2)) * t148) * t164;
	t1 = [(-t183 * t166 * t202 + (-t133 * t166 - t183 * t196) * t148) * t164, t113, 0, 0, 0; t153 * t126 * t205 + (-t135 * t126 + (t114 * t153 + t117 * t133) * t127) * t120 - ((t117 * t204 - t226 * t219) * t120 + (t117 * t205 + ((t118 * t148 * t199 + t202) * t217 + (0.2e1 * t199 * t221 - t118 + (-t184 * t164 + t118) * t148) * t216) * t120) * t127) * t183, (-t126 * t157 - t127 * t223) * t205 + (-t204 * t223 + t134 * t126 + (-t157 * t114 + t115 * t133 + (t169 * t206 - t113 * t153 - t119 * t135 + (t119 * t212 + t182) * t118) * t216 + (t118 * t119 * t153 - t136 + (t113 * t174 + (-qJD(2) * t119 - t118) * t171) * t169) * t217) * t127) * t120, 0, 0, 0; (t138 * t186 - t142 * t218) * t203 + ((t142 * qJD(3) - t136 * t170 + t173 * t198) * t138 - 0.2e1 * t142 * t201 + (t186 * t125 + (t186 * qJD(3) - t136 * t173 - t170 * t198) * t185 - t142 * t124) * t139) * t130, -t189 * t183 * t203 + (t189 * t133 - ((-qJD(3) * t138 + 0.2e1 * t201) * t173 + (t124 * t173 + (t125 + t227) * t170) * t139) * t183) * t130, -0.2e1 * t222 - 0.2e1 * (t124 * t139 * t130 - (-t130 * t220 - t139 * t222) * t185) * t185, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:42:21
	% EndTime: 2019-12-31 21:42:22
	% DurationCPUTime: 1.33s
	% Computational Cost: add. (4128->138), mult. (12381->282), div. (702->12), fcn. (15714->13), ass. (0->120)
	t224 = cos(pkin(5));
	t226 = sin(qJ(2));
	t292 = sin(qJ(1));
	t257 = t292 * t226;
	t247 = t224 * t257;
	t252 = qJD(2) * t292;
	t228 = cos(qJ(2));
	t229 = cos(qJ(1));
	t271 = t229 * t228;
	t222 = sin(pkin(5));
	t273 = t222 * t229;
	t297 = -qJD(1) * t247 - t226 * t252 + (qJD(2) * t224 + qJD(1)) * t271 - qJD(3) * t273;
	t210 = -t247 + t271;
	t225 = sin(qJ(3));
	t227 = cos(qJ(3));
	t258 = t222 * t292;
	t201 = t210 * t227 + t225 * t258;
	t221 = sin(pkin(10));
	t223 = cos(pkin(10));
	t256 = t292 * t228;
	t272 = t229 * t226;
	t238 = -t224 * t256 - t272;
	t177 = t201 * t221 + t223 * t238;
	t240 = -t224 * t272 - t256;
	t188 = t240 * qJD(1) + t238 * qJD(2);
	t239 = -t210 * t225 + t227 * t258;
	t255 = qJD(1) * t273;
	t167 = t239 * qJD(3) + t188 * t227 + t225 * t255;
	t251 = t292 * qJD(1);
	t259 = t224 * t271;
	t269 = qJD(2) * t228;
	t187 = -qJD(1) * t259 - t229 * t269 + (t224 * t252 + t251) * t226;
	t162 = t167 * t223 - t187 * t221;
	t178 = t201 * t223 - t221 * t238;
	t172 = 0.1e1 / t178;
	t173 = 0.1e1 / t178 ^ 2;
	t287 = t162 * t172 * t173;
	t250 = 0.2e1 * t177 * t287;
	t195 = -t225 * t240 + t227 * t273;
	t275 = t222 * t226;
	t205 = -t224 * t227 + t225 * t275;
	t184 = atan2(-t195, t205);
	t179 = sin(t184);
	t180 = cos(t184);
	t160 = -t179 * t195 + t180 * t205;
	t158 = 0.1e1 / t160 ^ 2;
	t192 = t239 ^ 2;
	t156 = t158 * t192 + 0.1e1;
	t166 = t201 * qJD(3) + t188 * t225 - t227 * t255;
	t286 = t166 * t158;
	t191 = t195 ^ 2;
	t203 = 0.1e1 / t205 ^ 2;
	t183 = t191 * t203 + 0.1e1;
	t181 = 0.1e1 / t183;
	t246 = t222 * t251;
	t268 = qJD(3) * t227;
	t168 = t297 * t225 - t227 * t246 - t240 * t268;
	t206 = t224 * t225 + t227 * t275;
	t254 = t222 * t269;
	t193 = t206 * qJD(3) + t225 * t254;
	t202 = 0.1e1 / t205;
	t278 = t195 * t203;
	t244 = -t168 * t202 + t193 * t278;
	t150 = t244 * t181;
	t245 = -t179 * t205 - t180 * t195;
	t145 = t245 * t150 - t168 * t179 + t180 * t193;
	t157 = 0.1e1 / t160;
	t159 = t157 * t158;
	t290 = t145 * t159;
	t267 = 0.2e1 * (-t192 * t290 - t239 * t286) / t156 ^ 2;
	t296 = t193 * t203;
	t207 = -t257 + t259;
	t274 = t222 * t228;
	t241 = -t202 * t207 + t274 * t278;
	t295 = t225 * t241;
	t169 = (qJD(3) * t240 + t246) * t225 + t297 * t227;
	t294 = -0.2e1 * t195;
	t293 = -0.2e1 * t239;
	t280 = t202 * t296;
	t289 = (t168 * t278 - t191 * t280) / t183 ^ 2;
	t288 = t158 * t239;
	t285 = t172 * t221;
	t284 = t173 * t177;
	t283 = t177 * t223;
	t282 = t179 * t239;
	t281 = t180 * t239;
	t279 = t195 * t202;
	t277 = t238 * t225;
	t276 = t238 * t227;
	t270 = qJD(2) * t226;
	t161 = t167 * t221 + t187 * t223;
	t171 = t177 ^ 2;
	t165 = t171 * t173 + 0.1e1;
	t266 = 0.2e1 * (t161 * t284 - t171 * t287) / t165 ^ 2;
	t265 = -0.2e1 * t289;
	t264 = t159 * t293;
	t263 = t202 * t289;
	t262 = t158 * t282;
	t261 = t158 * t281;
	t249 = t280 * t294;
	t197 = -t225 * t273 - t227 * t240;
	t243 = -t197 * t202 + t206 * t278;
	t242 = -qJD(3) * t277 + t187 * t227;
	t236 = -t179 + (t180 * t279 + t179) * t181;
	t194 = -t205 * qJD(3) + t227 * t254;
	t189 = t238 * qJD(1) + t240 * qJD(2);
	t186 = t210 * t221 + t223 * t276;
	t185 = -t210 * t223 + t221 * t276;
	t176 = -t197 * t223 + t207 * t221;
	t175 = -t197 * t221 - t207 * t223;
	t163 = 0.1e1 / t165;
	t154 = 0.1e1 / t156;
	t153 = t181 * t295;
	t152 = t243 * t181;
	t149 = t236 * t239;
	t147 = (-t179 * t207 + t180 * t274) * t225 + t245 * t153;
	t146 = t245 * t152 - t179 * t197 + t180 * t206;
	t144 = t243 * t265 + (t206 * t249 - t169 * t202 + (t168 * t206 + t193 * t197 + t194 * t195) * t203) * t181;
	t142 = t265 * t295 + (t241 * t268 + (t249 * t274 - t189 * t202 + (t193 * t207 + (t168 * t228 - t195 * t270) * t222) * t203) * t225) * t181;
	t1 = [t263 * t293 + (-t166 * t202 - t239 * t296) * t181, t142, t144, 0, 0; t195 * t157 * t267 + (-t168 * t157 + (t145 * t195 + t149 * t166) * t158) * t154 - (-t149 * t158 * t267 + (-0.2e1 * t149 * t290 + (-t150 * t181 * t279 + t265) * t262 + (t263 * t294 - t150 + (t150 - t244) * t181) * t261 - t236 * t286) * t154) * t239, (-t147 * t288 - t157 * t277) * t267 + (-t147 * t286 + (t187 * t225 + t238 * t268) * t157 + (t147 * t264 - t158 * t277) * t145 + (-t142 * t195 - t153 * t168 + (-t225 * t270 + t228 * t268) * t222 + (-t153 * t205 - t207 * t225) * t150) * t261 + (-t207 * t268 - t142 * t205 - t153 * t193 - t189 * t225 + (t153 * t195 - t225 * t274) * t150) * t262) * t154, (-t146 * t288 - t157 * t201) * t267 + (t146 * t145 * t264 + t167 * t157 + (-t201 * t145 - t146 * t166 + (-t144 * t195 - t152 * t168 + t194 + (-t152 * t205 - t197) * t150) * t281 + (-t144 * t205 - t152 * t193 - t169 + (t152 * t195 - t206) * t150) * t282) * t158) * t154, 0, 0; (-t172 * t175 + t176 * t284) * t266 + ((-t169 * t221 - t189 * t223) * t172 + t176 * t250 + (-t175 * t162 - (-t169 * t223 + t189 * t221) * t177 - t176 * t161) * t173) * t163, (-t172 * t185 + t186 * t284) * t266 + ((-t188 * t223 + t242 * t221) * t172 + t186 * t250 + (-t185 * t162 - (t188 * t221 + t242 * t223) * t177 - t186 * t161) * t173) * t163, -(-t173 * t283 + t285) * t239 * t266 + (t239 * t223 * t250 - t166 * t285 + (t166 * t283 - (t161 * t223 + t162 * t221) * t239) * t173) * t163, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 21:42:21
	% EndTime: 2019-12-31 21:42:22
	% DurationCPUTime: 1.41s
	% Computational Cost: add. (4943->149), mult. (13478->297), div. (726->12), fcn. (17045->13), ass. (0->127)
	t262 = cos(pkin(5));
	t264 = sin(qJ(2));
	t335 = sin(qJ(1));
	t297 = t335 * t264;
	t287 = t262 * t297;
	t291 = t335 * qJD(2);
	t266 = cos(qJ(2));
	t267 = cos(qJ(1));
	t313 = t267 * t266;
	t261 = sin(pkin(5));
	t315 = t261 * t267;
	t340 = -qJD(1) * t287 - t264 * t291 + (qJD(2) * t262 + qJD(1)) * t313 - qJD(3) * t315;
	t263 = sin(qJ(3));
	t265 = cos(qJ(3));
	t296 = t335 * t266;
	t314 = t267 * t264;
	t279 = -t262 * t314 - t296;
	t231 = -t263 * t279 + t265 * t315;
	t317 = t261 * t264;
	t242 = -t262 * t265 + t263 * t317;
	t220 = atan2(-t231, t242);
	t215 = sin(t220);
	t216 = cos(t220);
	t198 = -t215 * t231 + t216 * t242;
	t196 = 0.1e1 / t198 ^ 2;
	t247 = -t287 + t313;
	t298 = t261 * t335;
	t278 = -t247 * t263 + t265 * t298;
	t228 = t278 ^ 2;
	t194 = t228 * t196 + 0.1e1;
	t277 = -t262 * t296 - t314;
	t224 = qJD(1) * t279 + qJD(2) * t277;
	t237 = t247 * t265 + t263 * t298;
	t295 = qJD(1) * t315;
	t202 = qJD(3) * t237 + t224 * t263 - t265 * t295;
	t328 = t202 * t196;
	t227 = t231 ^ 2;
	t240 = 0.1e1 / t242 ^ 2;
	t219 = t227 * t240 + 0.1e1;
	t217 = 0.1e1 / t219;
	t292 = t335 * qJD(1);
	t286 = t261 * t292;
	t310 = qJD(3) * t265;
	t204 = t340 * t263 - t265 * t286 - t279 * t310;
	t243 = t262 * t263 + t265 * t317;
	t311 = qJD(2) * t266;
	t294 = t261 * t311;
	t229 = qJD(3) * t243 + t263 * t294;
	t239 = 0.1e1 / t242;
	t322 = t231 * t240;
	t283 = -t204 * t239 + t229 * t322;
	t186 = t283 * t217;
	t284 = -t215 * t242 - t216 * t231;
	t181 = t186 * t284 - t215 * t204 + t216 * t229;
	t195 = 0.1e1 / t198;
	t197 = t195 * t196;
	t333 = t181 * t197;
	t308 = 0.2e1 * (-t228 * t333 - t278 * t328) / t194 ^ 2;
	t339 = t229 * t240;
	t299 = t262 * t313;
	t244 = -t297 + t299;
	t316 = t261 * t266;
	t280 = -t239 * t244 + t316 * t322;
	t338 = t263 * t280;
	t205 = t263 * (qJD(3) * t279 + t286) + t340 * t265;
	t260 = pkin(10) + qJ(5);
	t258 = sin(t260);
	t259 = cos(t260);
	t214 = t237 * t259 - t258 * t277;
	t208 = 0.1e1 / t214;
	t209 = 0.1e1 / t214 ^ 2;
	t337 = -0.2e1 * t231;
	t336 = -0.2e1 * t278;
	t203 = qJD(3) * t278 + t224 * t265 + t263 * t295;
	t223 = -qJD(1) * t299 - t267 * t311 + (t262 * t291 + t292) * t264;
	t189 = qJD(5) * t214 + t203 * t258 + t223 * t259;
	t213 = t237 * t258 + t259 * t277;
	t207 = t213 ^ 2;
	t201 = t207 * t209 + 0.1e1;
	t327 = t209 * t213;
	t309 = qJD(5) * t213;
	t190 = t203 * t259 - t223 * t258 - t309;
	t330 = t190 * t208 * t209;
	t332 = (t189 * t327 - t207 * t330) / t201 ^ 2;
	t324 = t239 * t339;
	t331 = (t204 * t322 - t227 * t324) / t219 ^ 2;
	t329 = t196 * t278;
	t326 = t215 * t278;
	t325 = t216 * t278;
	t323 = t231 * t239;
	t321 = t277 * t263;
	t320 = t277 * t265;
	t319 = t258 * t208;
	t318 = t259 * t213;
	t312 = qJD(2) * t264;
	t307 = -0.2e1 * t332;
	t306 = 0.2e1 * t332;
	t305 = -0.2e1 * t331;
	t304 = t197 * t336;
	t303 = t239 * t331;
	t302 = t196 * t326;
	t301 = t196 * t325;
	t300 = t213 * t330;
	t290 = 0.2e1 * t300;
	t289 = t324 * t337;
	t233 = -t263 * t315 - t265 * t279;
	t285 = -qJD(5) * t320 + t224;
	t212 = -t233 * t259 + t244 * t258;
	t211 = -t233 * t258 - t244 * t259;
	t282 = t209 * t318 - t319;
	t281 = -t233 * t239 + t243 * t322;
	t275 = -t215 + (t216 * t323 + t215) * t217;
	t274 = -qJD(3) * t321 + qJD(5) * t247 + t223 * t265;
	t230 = -qJD(3) * t242 + t265 * t294;
	t225 = qJD(1) * t277 + qJD(2) * t279;
	t222 = t247 * t258 + t259 * t320;
	t221 = -t247 * t259 + t258 * t320;
	t199 = 0.1e1 / t201;
	t192 = 0.1e1 / t194;
	t191 = t217 * t338;
	t188 = t281 * t217;
	t185 = t275 * t278;
	t183 = (-t215 * t244 + t216 * t316) * t263 + t284 * t191;
	t182 = t188 * t284 - t215 * t233 + t216 * t243;
	t180 = t281 * t305 + (t243 * t289 - t205 * t239 + (t204 * t243 + t229 * t233 + t230 * t231) * t240) * t217;
	t178 = t305 * t338 + (t280 * t310 + (t289 * t316 - t225 * t239 + (t229 * t244 + (t204 * t266 - t231 * t312) * t261) * t240) * t263) * t217;
	t1 = [t303 * t336 + (-t202 * t239 - t278 * t339) * t217, t178, t180, 0, 0; t231 * t195 * t308 + (-t204 * t195 + (t181 * t231 + t185 * t202) * t196) * t192 - (-t185 * t196 * t308 + (-0.2e1 * t185 * t333 + (-t186 * t217 * t323 + t305) * t302 + (t303 * t337 - t186 + (t186 - t283) * t217) * t301 - t275 * t328) * t192) * t278, (-t183 * t329 - t195 * t321) * t308 + (-t183 * t328 + (t223 * t263 + t277 * t310) * t195 + (t183 * t304 - t196 * t321) * t181 + (-t178 * t231 - t191 * t204 + (-t263 * t312 + t266 * t310) * t261 + (-t191 * t242 - t244 * t263) * t186) * t301 + (-t244 * t310 - t178 * t242 - t191 * t229 - t225 * t263 + (t191 * t231 - t263 * t316) * t186) * t302) * t192, (-t182 * t329 - t195 * t237) * t308 + (t182 * t181 * t304 + t203 * t195 + (-t237 * t181 - t182 * t202 + (-t180 * t231 - t188 * t204 + t230 + (-t188 * t242 - t233) * t186) * t325 + (-t180 * t242 - t188 * t229 - t205 + (t188 * t231 - t243) * t186) * t326) * t196) * t192, 0, 0; (-t208 * t211 + t212 * t327) * t306 + ((qJD(5) * t212 - t205 * t258 - t225 * t259) * t208 + t212 * t290 + (-t211 * t190 - (-qJD(5) * t211 - t205 * t259 + t225 * t258) * t213 - t212 * t189) * t209) * t199, (-t208 * t221 + t222 * t327) * t306 + (t222 * t290 - t285 * t208 * t259 + t274 * t319 + (-t213 * t258 * t285 - t222 * t189 - t221 * t190 - t274 * t318) * t209) * t199, -t282 * t278 * t307 + (t282 * t202 - ((-qJD(5) * t208 - 0.2e1 * t300) * t259 + (t189 * t259 + (t190 - t309) * t258) * t209) * t278) * t199, 0, t307 + 0.2e1 * (t189 * t209 * t199 + (-t199 * t330 - t209 * t332) * t213) * t213;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end