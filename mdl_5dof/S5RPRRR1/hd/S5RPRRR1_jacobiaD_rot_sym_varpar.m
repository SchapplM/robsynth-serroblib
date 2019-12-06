% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRR1
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
%   Wie in S5RPRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiaD_rot_sym_varpar: pkin has to be [1x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:18
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:19
	% DurationCPUTime: 0.82s
	% Computational Cost: add. (1002->94), mult. (2519->208), div. (480->12), fcn. (2968->9), ass. (0->94)
	t121 = sin(qJ(4));
	t123 = sin(qJ(1));
	t124 = cos(qJ(4));
	t125 = cos(qJ(3));
	t126 = cos(qJ(1));
	t163 = t125 * t126;
	t103 = t121 * t163 - t123 * t124;
	t192 = 0.2e1 * t103;
	t116 = t123 ^ 2;
	t122 = sin(qJ(3));
	t115 = t122 ^ 2;
	t118 = 0.1e1 / t125 ^ 2;
	t170 = t115 * t118;
	t110 = t116 * t170 + 0.1e1;
	t114 = t122 * t115;
	t117 = 0.1e1 / t125;
	t169 = t117 * t122;
	t134 = qJD(3) * (t114 * t117 * t118 + t169);
	t161 = qJD(1) * t126;
	t147 = t123 * t161;
	t175 = 0.1e1 / t110 ^ 2 * (t116 * t134 + t147 * t170);
	t191 = -0.2e1 * t175;
	t108 = 0.1e1 / t110;
	t142 = 0.1e1 + t170;
	t187 = t123 * t142;
	t95 = t108 * t187;
	t190 = t123 * t95 - 0.1e1;
	t140 = -qJD(1) * t125 + qJD(4);
	t157 = qJD(3) * t126;
	t145 = t122 * t157;
	t168 = t121 * t126;
	t156 = qJD(4) * t125;
	t186 = qJD(1) - t156;
	t84 = t186 * t168 + (t140 * t123 - t145) * t124;
	t166 = t123 * t121;
	t104 = t124 * t163 + t166;
	t99 = 0.1e1 / t104 ^ 2;
	t189 = t84 * t99;
	t177 = t103 * t99;
	t98 = 0.1e1 / t104;
	t136 = -t121 * t98 + t124 * t177;
	t97 = t103 ^ 2;
	t96 = t97 * t99 + 0.1e1;
	t93 = 0.1e1 / t96;
	t188 = t136 * t93;
	t165 = t123 * t122;
	t107 = atan2(-t165, -t125);
	t106 = cos(t107);
	t105 = sin(t107);
	t150 = t105 * t165;
	t91 = -t106 * t125 - t150;
	t88 = 0.1e1 / t91;
	t89 = 0.1e1 / t91 ^ 2;
	t185 = t108 - 0.1e1;
	t120 = t126 ^ 2;
	t158 = qJD(3) * t125;
	t151 = t89 * t158;
	t148 = t122 * t161;
	t159 = qJD(3) * t123;
	t171 = t106 * t122;
	t146 = t118 * t159;
	t82 = (-(-t123 * t158 - t148) * t117 + t115 * t146) * t108;
	t77 = (-t123 * t82 + qJD(3)) * t171 + (-t148 + (t82 - t159) * t125) * t105;
	t183 = t77 * t88 * t89;
	t87 = t115 * t120 * t89 + 0.1e1;
	t184 = (t120 * t122 * t151 + (-t120 * t183 - t89 * t147) * t115) / t87 ^ 2;
	t178 = t98 * t189;
	t164 = t123 * t125;
	t135 = t121 * t164 + t124 * t126;
	t144 = t124 * t156;
	t83 = t135 * qJD(1) - qJD(4) * t166 + t121 * t145 - t126 * t144;
	t182 = (-t83 * t177 - t97 * t178) / t96 ^ 2;
	t181 = t83 * t99;
	t85 = 0.1e1 / t87;
	t180 = t85 * t89;
	t179 = t88 * t85;
	t172 = qJD(3) * t95;
	t167 = t122 * t126;
	t162 = qJD(1) * t123;
	t160 = qJD(3) * t122;
	t155 = 0.2e1 * t183;
	t154 = -0.2e1 * t182;
	t153 = t88 * t184;
	t149 = t108 * t115 * t117;
	t143 = 0.2e1 * t89 * t184;
	t141 = t117 * t191;
	t139 = t123 * t149;
	t138 = t142 * t126;
	t137 = t178 * t192 + t181;
	t102 = -t124 * t164 + t168;
	t81 = (t185 * t122 * t105 - t106 * t139) * t126;
	t80 = -t190 * t171 + (-t123 + t95) * t125 * t105;
	t78 = t187 * t191 + (qJD(1) * t138 + 0.2e1 * t123 * t134) * t108;
	t1 = [t141 * t167 + (qJD(3) * t138 - t162 * t169) * t108, 0, t78, 0, 0; (-t158 * t179 + (0.2e1 * t153 + (qJD(1) * t81 + t77) * t180) * t122) * t123 + (t81 * t143 * t122 + (-t81 * t151 + (t81 * t155 + ((0.2e1 * t122 * t175 - t82 * t139 - t185 * t158) * t105 + (t115 * t123 * t141 + t122 * t82 + (t114 * t146 - (t82 - 0.2e1 * t159) * t122) * t108) * t106) * t89 * t126) * t122 + (-t88 + (-(t116 - t120) * t106 * t149 + t185 * t150) * t89) * t122 * qJD(1)) * t85) * t126, 0, (-t162 * t179 + (-0.2e1 * t153 + (-qJD(3) * t80 - t77) * t180) * t126) * t125 + (t80 * t126 * t143 + (-t88 * t157 - ((-t123 * t78 - t161 * t95) * t106 + (t190 * t82 + t159 - t172) * t105) * t89 * t167 + (t126 * t155 + t89 * t162) * t80) * t85 - ((t78 - t161) * t105 + (t82 * t95 + qJD(3) + (-t82 - t172) * t123) * t106) * t163 * t180) * t122, 0, 0; 0.2e1 * (t102 * t177 + t135 * t98) * t182 + (t135 * t189 + t137 * t102 + ((qJD(1) * t124 + t121 * t160 - t144) * t98 - (-t186 * t121 + t124 * t160) * t177) * t123 - t136 * t126 * t140) * t93, 0, t125 * t157 * t188 + (-t162 * t188 + (t136 * t154 + ((-qJD(4) * t103 + t84) * t99 * t121 + (-qJD(4) * t98 - t137) * t124) * t93) * t126) * t122, t154 + (-t93 * t181 + (-t93 * t178 - t99 * t182) * t103) * t192, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:10:18
	% EndTime: 2019-12-05 18:10:20
	% DurationCPUTime: 1.48s
	% Computational Cost: add. (2288->151), mult. (7495->320), div. (1129->14), fcn. (9361->11), ass. (0->133)
	t200 = cos(qJ(3));
	t196 = sin(qJ(4));
	t287 = sin(qJ(1));
	t237 = t287 * t196;
	t199 = cos(qJ(4));
	t201 = cos(qJ(1));
	t263 = t201 * t199;
	t177 = t200 * t263 + t237;
	t195 = sin(qJ(5));
	t198 = cos(qJ(5));
	t197 = sin(qJ(3));
	t265 = t197 * t201;
	t161 = t177 * t195 - t198 * t265;
	t291 = 0.2e1 * t161;
	t162 = t177 * t198 + t195 * t265;
	t156 = 0.1e1 / t162;
	t157 = 0.1e1 / t162 ^ 2;
	t278 = t157 * t161;
	t217 = t156 * t195 - t198 * t278;
	t290 = -qJD(5) * t199 + qJD(3);
	t173 = t200 * t237 + t263;
	t228 = t287 * qJD(4);
	t220 = t196 * t228;
	t257 = qJD(4) * t201;
	t232 = t199 * t257;
	t260 = qJD(3) * t201;
	t234 = t197 * t260;
	t151 = t173 * qJD(1) + t196 * t234 - t200 * t232 - t220;
	t236 = t287 * t199;
	t264 = t201 * t196;
	t176 = t200 * t264 - t236;
	t189 = 0.1e1 / t196;
	t190 = 0.1e1 / t196 ^ 2;
	t192 = 0.1e1 / t197;
	t193 = 0.1e1 / t197 ^ 2;
	t261 = qJD(3) * t200;
	t235 = t193 * t261;
	t258 = qJD(4) * t199;
	t270 = t189 * t192;
	t289 = (t190 * t192 * t258 + t189 * t235) * t176 + t151 * t270;
	t267 = t197 * t196;
	t167 = atan2(-t173, t267);
	t164 = cos(t167);
	t163 = sin(t167);
	t276 = t163 * t173;
	t147 = t164 * t267 - t276;
	t144 = 0.1e1 / t147;
	t145 = 0.1e1 / t147 ^ 2;
	t288 = 0.2e1 * t176;
	t171 = t173 ^ 2;
	t269 = t190 * t193;
	t168 = t171 * t269 + 0.1e1;
	t165 = 0.1e1 / t168;
	t214 = t196 * t261 + t197 * t258;
	t241 = t173 * t269;
	t229 = qJD(3) * t287;
	t221 = t197 * t229;
	t230 = qJD(1) * t287;
	t262 = qJD(1) * t201;
	t153 = (t228 * t200 - t230) * t199 + (t262 * t200 - t221 - t257) * t196;
	t244 = t153 * t270;
	t135 = (t214 * t241 - t244) * t165;
	t212 = -t135 * t173 + t214;
	t130 = (-t135 * t267 - t153) * t163 + t212 * t164;
	t146 = t144 * t145;
	t286 = t130 * t146;
	t211 = qJD(5) * t177 + t197 * t230 - t200 * t260;
	t152 = (-qJD(4) * t200 + qJD(1)) * t264 + (-t200 * t230 + t228 - t234) * t199;
	t219 = qJD(5) * t265 + t152;
	t137 = t219 * t195 + t211 * t198;
	t155 = t161 ^ 2;
	t150 = t155 * t157 + 0.1e1;
	t138 = -t211 * t195 + t219 * t198;
	t158 = t156 * t157;
	t283 = t138 * t158;
	t285 = (t137 * t278 - t155 * t283) / t150 ^ 2;
	t191 = t189 * t190;
	t194 = t192 * t193;
	t233 = t193 * t258;
	t284 = (t153 * t241 + (-t190 * t194 * t261 - t191 * t233) * t171) / t168 ^ 2;
	t282 = t145 * t151;
	t281 = t145 * t176;
	t148 = 0.1e1 / t150;
	t280 = t148 * t157;
	t266 = t197 * t199;
	t170 = (t195 * t200 - t198 * t266) * t201;
	t277 = t161 * t170;
	t275 = t163 * t176;
	t274 = t163 * t197;
	t273 = t164 * t173;
	t272 = t164 * t176;
	t271 = t164 * t200;
	t268 = t190 * t199;
	t259 = qJD(4) * t196;
	t172 = t176 ^ 2;
	t142 = t145 * t172 + 0.1e1;
	t255 = 0.2e1 * (-t151 * t281 - t172 * t286) / t142 ^ 2;
	t254 = -0.2e1 * t285;
	t253 = 0.2e1 * t285;
	t252 = -0.2e1 * t284;
	t251 = t146 * t288;
	t250 = t157 * t285;
	t249 = t192 * t284;
	t248 = t137 * t280;
	t247 = t161 * t283;
	t246 = t145 * t275;
	t242 = t173 * t270;
	t240 = t189 * t193 * t200;
	t239 = t197 * t287;
	t238 = t200 * t287;
	t215 = t173 * t240 + t287;
	t143 = t215 * t165;
	t231 = t287 - t143;
	t227 = t144 * t255;
	t226 = t145 * t255;
	t224 = t189 * t249;
	t223 = t195 * t239;
	t222 = t198 * t239;
	t154 = t177 * qJD(1) - t199 * t221 - t200 * t220 - t232;
	t218 = -qJD(5) * t239 - t154;
	t175 = t200 * t236 - t264;
	t216 = t173 * t268 - t175 * t189;
	t210 = -qJD(5) * t175 + t197 * t262 + t200 * t229;
	t169 = (-t195 * t266 - t198 * t200) * t201;
	t160 = -t175 * t198 - t223;
	t140 = 0.1e1 / t142;
	t139 = t216 * t192 * t165;
	t134 = (-t163 + (t164 * t242 + t163) * t165) * t176;
	t133 = -t143 * t273 + (t231 * t274 + t271) * t196;
	t131 = t164 * t266 - t163 * t175 + (-t163 * t267 - t273) * t139;
	t129 = t215 * t252 + (t153 * t240 + t262 + (-t190 * t200 * t233 + (-0.2e1 * t194 * t200 ^ 2 - t192) * t189 * qJD(3)) * t173) * t165;
	t127 = -0.2e1 * t216 * t249 + (-t216 * t235 + (t153 * t268 - t154 * t189 + (t175 * t268 + (-0.2e1 * t191 * t199 ^ 2 - t189) * t173) * qJD(4)) * t192) * t165;
	t1 = [t289 * t165 + t224 * t288, 0, t129, t127, 0; t173 * t227 + (-t153 * t144 + (t130 * t173 + t134 * t151) * t145) * t140 + (t134 * t226 + (0.2e1 * t134 * t286 + (t151 * t165 - t151 - (-t135 * t165 * t242 + t252) * t176) * t145 * t163 + (-(-0.2e1 * t173 * t224 - t135) * t281 + (-(t135 + t244) * t176 + t289 * t173) * t145 * t165) * t164) * t140) * t176, 0, t133 * t176 * t226 + (-(-t129 * t273 + (t135 * t276 - t153 * t164) * t143) * t281 + (t130 * t251 + t282) * t133 + (-t144 * t265 - (-t143 * t274 + t163 * t239 + t271) * t281) * t258) * t140 + (t227 * t265 + ((-t144 * t260 - (t231 * qJD(3) - t135) * t246) * t200 + (t144 * t230 + (t201 * t130 - (-t129 + t262) * t275 - (t231 * t135 - qJD(3)) * t272) * t145) * t197) * t140) * t196, (t131 * t281 - t144 * t177) * t255 + (t131 * t282 + t152 * t144 + (t131 * t251 - t145 * t177) * t130 - (t199 * t261 - t197 * t259 - t127 * t173 - t139 * t153 + (-t139 * t267 - t175) * t135) * t145 * t272 - (-t154 + (-t127 * t196 - t135 * t199) * t197 - t212 * t139) * t246) * t140, 0; (t250 * t291 - t248) * t160 + (-t138 * t280 + t156 * t254) * (-t175 * t195 + t222) + ((t218 * t195 + t210 * t198) * t156 - (-t210 * t195 + t218 * t198) * t278 + 0.2e1 * t160 * t247) * t148, 0, (-t156 * t169 + t157 * t277) * t253 + (-t170 * t137 * t157 + (-t157 * t169 + 0.2e1 * t158 * t277) * t138 + ((t198 * t238 + t199 * t223) * t156 - (-t195 * t238 + t199 * t222) * t278) * qJD(1) + (((t195 * t259 + t290 * t198) * t156 - (-t290 * t195 + t198 * t259) * t278) * t197 + t217 * t200 * (-qJD(3) * t199 + qJD(5))) * t201) * t148, t217 * t176 * t253 + (t217 * t151 + ((-qJD(5) * t156 - 0.2e1 * t247) * t198 + (t137 * t198 + (-qJD(5) * t161 + t138) * t195) * t157) * t176) * t148, t254 + (t248 + (-t148 * t283 - t250) * t161) * t291;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end