% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6PRPPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:24
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiaD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:27
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (46->7), mult. (159->21), div. (18->4), fcn. (175->5), ass. (0->15)
	t39 = cos(pkin(10));
	t41 = sin(qJ(2));
	t42 = cos(qJ(2));
	t45 = sin(pkin(10)) * cos(pkin(6));
	t37 = t39 * t42 - t41 * t45;
	t34 = 0.1e1 / t37 ^ 2;
	t49 = qJD(2) * t34;
	t36 = t39 * t41 + t42 * t45;
	t33 = t36 ^ 2;
	t30 = t33 * t34 + 0.1e1;
	t46 = t37 * t49;
	t47 = t36 / t37 * t49;
	t48 = (t33 * t47 + t36 * t46) / t30 ^ 2;
	t28 = 0.1e1 / t30;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, -0.2e1 * t48 + 0.2e1 * (t28 * t46 + (t28 * t47 - t34 * t48) * t36) * t36, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:27
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (153->13), mult. (451->28), div. (25->4), fcn. (552->7), ass. (0->20)
	t63 = sin(pkin(11));
	t64 = cos(pkin(11));
	t65 = cos(pkin(10));
	t66 = sin(qJ(2));
	t67 = cos(qJ(2));
	t78 = cos(pkin(6));
	t74 = t67 * t78;
	t75 = t66 * t78;
	t77 = sin(pkin(10));
	t54 = (t63 * t75 - t64 * t74) * t77 + t65 * (-t67 * t63 - t66 * t64);
	t48 = t54 * qJD(2);
	t72 = (-t63 * t74 - t64 * t75) * t77 - t65 * (t66 * t63 - t67 * t64);
	t51 = 0.1e1 / t72 ^ 2;
	t84 = t51 * t54 ^ 2;
	t50 = 0.1e1 / t72;
	t83 = t50 * t84;
	t82 = t51 * t72;
	t76 = t82 * t48;
	t46 = 0.1e1 + t84;
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0.2e1 * (-t50 * t72 - t84) / t46 ^ 2 * (-t48 * t83 - t76) + (-0.2e1 * t76 + (t50 - t82 - 0.2e1 * t83) * t48) / t46, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:28
	% DurationCPUTime: 0.59s
	% Computational Cost: add. (1475->50), mult. (4551->125), div. (254->12), fcn. (5909->13), ass. (0->70)
	t150 = cos(pkin(6));
	t144 = sin(pkin(11));
	t148 = cos(pkin(11));
	t151 = sin(qJ(2));
	t152 = cos(qJ(2));
	t162 = t152 * t144 + t151 * t148;
	t136 = t162 * t150;
	t139 = t151 * t144 - t152 * t148;
	t137 = t139 * qJD(2);
	t145 = sin(pkin(10));
	t149 = cos(pkin(10));
	t160 = t139 * t150;
	t121 = -t145 * t162 - t149 * t160;
	t146 = sin(pkin(6));
	t134 = t139 * t146;
	t107 = atan2(t121, t134);
	t102 = sin(t107);
	t103 = cos(t107);
	t98 = t102 * t121 + t103 * t134;
	t95 = 0.1e1 / t98;
	t143 = sin(pkin(12));
	t147 = cos(pkin(12));
	t163 = -t145 * t136 - t149 * t139;
	t168 = t145 * t146;
	t113 = t143 * t168 + t147 * t163;
	t109 = 0.1e1 / t113;
	t131 = 0.1e1 / t134;
	t110 = 0.1e1 / t113 ^ 2;
	t132 = 0.1e1 / t134 ^ 2;
	t96 = 0.1e1 / t98 ^ 2;
	t112 = t143 * t163 - t147 * t168;
	t108 = t112 ^ 2;
	t101 = t108 * t110 + 0.1e1;
	t111 = t109 * t110;
	t130 = t150 * t137;
	t138 = t162 * qJD(2);
	t117 = t145 * t130 - t149 * t138;
	t176 = 0.1e1 / t101 ^ 2 * (-t108 * t111 * t147 + t110 * t112 * t143) * t117;
	t99 = 0.1e1 / t101;
	t175 = t117 * t99;
	t123 = t145 * t160 - t149 * t162;
	t174 = t123 * t96;
	t173 = t109 * t143;
	t172 = t112 * t147;
	t171 = t121 * t132;
	t135 = t162 * t146;
	t170 = t121 * t135;
	t128 = qJD(2) * t135;
	t169 = t128 * t131 * t132;
	t120 = -t149 * t136 + t145 * t139;
	t118 = t121 ^ 2;
	t106 = t118 * t132 + 0.1e1;
	t104 = 0.1e1 / t106;
	t161 = -t120 * t131 + t132 * t170;
	t90 = t161 * t104;
	t165 = t134 * t90 + t120;
	t164 = -t121 * t90 + t135;
	t159 = qJD(2) * t136;
	t129 = t146 * t137;
	t119 = t123 ^ 2;
	t116 = t149 * t137 + t145 * t159;
	t115 = t149 * t130 + t145 * t138;
	t114 = t145 * t137 - t149 * t159;
	t97 = t95 * t96;
	t93 = t119 * t96 + 0.1e1;
	t89 = (t114 * t131 - t128 * t171) * t104;
	t87 = t165 * t102 + t164 * t103;
	t86 = (t121 * t89 + t128) * t103 + (-t134 * t89 + t114) * t102;
	t85 = 0.2e1 * t161 * (t114 * t171 - t118 * t169) / t106 ^ 2 + (0.2e1 * t169 * t170 + t115 * t131 + (-t114 * t135 - t120 * t128 + t121 * t129) * t132) * t104;
	t1 = [0, t85, 0, 0, 0, 0; 0, 0.2e1 * (-t163 * t95 - t87 * t174) / t93 ^ 2 * (-t119 * t97 * t86 + t116 * t174) + (t117 * t95 + t87 * t116 * t96 + (-0.2e1 * t87 * t123 * t97 - t163 * t96) * t86 + ((-t114 * t90 + t121 * t85 + t165 * t89 - t129) * t103 + (t128 * t90 - t134 * t85 - t164 * t89 + t115) * t102) * t174) / t93, 0, 0, 0, 0; 0, (-t110 * t172 + t173) * t99 * t116 + 0.2e1 * (-t173 * t176 + (t111 * t172 * t175 + (t112 * t176 - t143 * t175) * t110) * t147) * t123, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:29
	% DurationCPUTime: 0.73s
	% Computational Cost: add. (1932->59), mult. (5333->135), div. (281->12), fcn. (6885->13), ass. (0->72)
	t188 = pkin(12) + qJ(5);
	t186 = sin(t188);
	t187 = cos(t188);
	t194 = cos(pkin(6));
	t189 = sin(pkin(11));
	t192 = cos(pkin(11));
	t195 = sin(qJ(2));
	t196 = cos(qJ(2));
	t208 = t189 * t196 + t195 * t192;
	t178 = t208 * t194;
	t181 = t195 * t189 - t192 * t196;
	t190 = sin(pkin(10));
	t193 = cos(pkin(10));
	t209 = -t178 * t190 - t181 * t193;
	t191 = sin(pkin(6));
	t212 = t190 * t191;
	t205 = -t186 * t209 + t187 * t212;
	t222 = t205 * qJD(5);
	t179 = t181 * qJD(2);
	t204 = t181 * t194;
	t163 = -t190 * t208 - t193 * t204;
	t176 = t181 * t191;
	t149 = atan2(t163, t176);
	t144 = sin(t149);
	t145 = cos(t149);
	t138 = t144 * t163 + t145 * t176;
	t135 = 0.1e1 / t138;
	t155 = t186 * t212 + t187 * t209;
	t151 = 0.1e1 / t155;
	t173 = 0.1e1 / t176;
	t136 = 0.1e1 / t138 ^ 2;
	t152 = 0.1e1 / t155 ^ 2;
	t174 = 0.1e1 / t176 ^ 2;
	t160 = t163 ^ 2;
	t148 = t160 * t174 + 0.1e1;
	t146 = 0.1e1 / t148;
	t203 = qJD(2) * t178;
	t156 = t190 * t179 - t193 * t203;
	t177 = t208 * t191;
	t170 = qJD(2) * t177;
	t216 = t163 * t174;
	t129 = (t156 * t173 - t170 * t216) * t146;
	t210 = -t144 * t176 + t145 * t163;
	t126 = t210 * t129 + t144 * t156 + t145 * t170;
	t221 = t126 * t135 * t136;
	t150 = t205 ^ 2;
	t141 = t150 * t152 + 0.1e1;
	t172 = t194 * t179;
	t180 = t208 * qJD(2);
	t159 = t172 * t190 - t180 * t193;
	t142 = t155 * qJD(5) + t159 * t186;
	t217 = t152 * t205;
	t143 = t159 * t187 + t222;
	t218 = t143 * t151 * t152;
	t220 = (-t142 * t217 - t150 * t218) / t141 ^ 2;
	t165 = t190 * t204 - t193 * t208;
	t219 = t136 * t165;
	t215 = t163 * t177;
	t214 = t170 * t173 * t174;
	t207 = -t151 * t186 - t187 * t217;
	t162 = -t178 * t193 + t181 * t190;
	t206 = -t162 * t173 + t174 * t215;
	t171 = t191 * t179;
	t161 = t165 ^ 2;
	t158 = t179 * t193 + t190 * t203;
	t157 = t172 * t193 + t180 * t190;
	t139 = 0.1e1 / t141;
	t133 = t136 * t161 + 0.1e1;
	t130 = t206 * t146;
	t127 = -t210 * t130 + t144 * t162 + t145 * t177;
	t125 = 0.2e1 * t206 / t148 ^ 2 * (t156 * t216 - t160 * t214) + (0.2e1 * t214 * t215 + t157 * t173 + (-t156 * t177 - t162 * t170 + t163 * t171) * t174) * t146;
	t1 = [0, t125, 0, 0, 0, 0; 0, 0.2e1 * (-t127 * t219 - t135 * t209) / t133 ^ 2 * (t158 * t219 - t161 * t221) + (t159 * t135 + (-t126 * t209 + t127 * t158) * t136 + (-0.2e1 * t127 * t221 + ((t125 * t163 - t130 * t156 - t171 + (t130 * t176 + t162) * t129) * t145 + (-t125 * t176 + t130 * t170 + t157 + (t130 * t163 - t177) * t129) * t144) * t136) * t165) / t133, 0, 0, 0, 0; 0, 0.2e1 * t207 * t165 * t220 + (-t207 * t158 + ((qJD(5) * t151 - 0.2e1 * t205 * t218) * t187 + (-t142 * t187 + (-t143 - t222) * t186) * t152) * t165) * t139, 0, 0, -0.2e1 * t220 - 0.2e1 * (t139 * t142 * t152 - (-t139 * t218 - t152 * t220) * t205) * t205, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:24:28
	% EndTime: 2019-10-09 21:24:30
	% DurationCPUTime: 1.85s
	% Computational Cost: add. (8444->115), mult. (16325->235), div. (559->12), fcn. (21250->15), ass. (0->112)
	t281 = sin(pkin(11));
	t284 = cos(pkin(11));
	t288 = sin(qJ(2));
	t290 = cos(qJ(2));
	t271 = t288 * t281 - t290 * t284;
	t286 = cos(pkin(6));
	t299 = t271 * t286;
	t265 = qJD(2) * t299;
	t304 = t290 * t281 + t288 * t284;
	t270 = t304 * qJD(2);
	t282 = sin(pkin(10));
	t285 = cos(pkin(10));
	t246 = -t285 * t265 - t282 * t270;
	t268 = t304 * t286;
	t252 = t285 * t268 - t282 * t271;
	t280 = pkin(12) + qJ(5);
	t278 = sin(t280);
	t283 = sin(pkin(6));
	t320 = t283 * t285;
	t310 = t278 * t320;
	t279 = cos(t280);
	t316 = qJD(5) * t279;
	t218 = -qJD(5) * t310 + t246 * t278 + t252 * t316;
	t240 = t252 * t278 + t279 * t320;
	t238 = t240 ^ 2;
	t267 = t304 * t283;
	t259 = t267 * t278 - t286 * t279;
	t257 = 0.1e1 / t259 ^ 2;
	t232 = t238 * t257 + 0.1e1;
	t230 = 0.1e1 / t232;
	t260 = t267 * t279 + t286 * t278;
	t266 = t271 * t283;
	t264 = qJD(2) * t266;
	t236 = qJD(5) * t260 - t264 * t278;
	t256 = 0.1e1 / t259;
	t325 = t240 * t257;
	t202 = (-t218 * t256 + t236 * t325) * t230;
	t233 = atan2(-t240, t259);
	t228 = sin(t233);
	t229 = cos(t233);
	t307 = -t228 * t259 - t229 * t240;
	t198 = t202 * t307 - t228 * t218 + t229 * t236;
	t212 = -t228 * t240 + t229 * t259;
	t209 = 0.1e1 / t212;
	t210 = 0.1e1 / t212 ^ 2;
	t339 = t198 * t209 * t210;
	t305 = -t282 * t268 - t285 * t271;
	t321 = t282 * t283;
	t300 = -t278 * t305 + t279 * t321;
	t338 = -0.2e1 * t300 * t339;
	t251 = -t282 * t304 - t285 * t299;
	t301 = -t251 * t256 - t266 * t325;
	t337 = t278 * t301;
	t326 = t236 * t256 * t257;
	t336 = -0.2e1 * (t218 * t325 - t238 * t326) / t232 ^ 2;
	t244 = t278 * t321 + t279 * t305;
	t289 = cos(qJ(6));
	t254 = t282 * t299 - t285 * t304;
	t287 = sin(qJ(6));
	t323 = t254 * t287;
	t227 = t244 * t289 - t323;
	t223 = 0.1e1 / t227;
	t224 = 0.1e1 / t227 ^ 2;
	t306 = t282 * t265 - t285 * t270;
	t221 = qJD(5) * t300 + t279 * t306;
	t269 = t271 * qJD(2);
	t298 = t286 * t270;
	t247 = t285 * t269 + t282 * t298;
	t213 = qJD(6) * t227 + t221 * t287 + t247 * t289;
	t322 = t254 * t289;
	t226 = t244 * t287 + t322;
	t222 = t226 ^ 2;
	t217 = t222 * t224 + 0.1e1;
	t330 = t224 * t226;
	t315 = qJD(6) * t226;
	t214 = t221 * t289 - t247 * t287 - t315;
	t333 = t214 * t223 * t224;
	t335 = (t213 * t330 - t222 * t333) / t217 ^ 2;
	t334 = t210 * t300;
	t220 = qJD(5) * t244 + t278 * t306;
	t332 = t220 * t210;
	t331 = t223 * t287;
	t329 = t226 * t289;
	t328 = t228 * t300;
	t327 = t229 * t300;
	t324 = t254 * t278;
	t239 = t300 ^ 2;
	t208 = t239 * t210 + 0.1e1;
	t314 = 0.2e1 * (-t239 * t339 - t300 * t332) / t208 ^ 2;
	t313 = -0.2e1 * t335;
	t311 = t226 * t333;
	t309 = -0.2e1 * t240 * t326;
	t308 = qJD(6) * t254 * t279 - t306;
	t303 = t329 * t224 - t331;
	t242 = t252 * t279 - t310;
	t302 = -t242 * t256 + t260 * t325;
	t297 = -qJD(5) * t324 + qJD(6) * t305 + t247 * t279;
	t263 = t283 * t270;
	t245 = t282 * t269 - t285 * t298;
	t237 = -qJD(5) * t259 - t264 * t279;
	t235 = t279 * t322 + t287 * t305;
	t234 = t279 * t323 - t289 * t305;
	t219 = -qJD(5) * t240 + t246 * t279;
	t215 = 0.1e1 / t217;
	t206 = 0.1e1 / t208;
	t204 = t230 * t337;
	t203 = t302 * t230;
	t200 = (-t228 * t251 - t229 * t266) * t278 + t307 * t204;
	t199 = t203 * t307 - t228 * t242 + t229 * t260;
	t196 = t302 * t336 + (t260 * t309 - t219 * t256 + (t218 * t260 + t236 * t242 + t237 * t240) * t257) * t230;
	t195 = t336 * t337 + (t301 * t316 + (-t266 * t309 - t245 * t256 + (-t218 * t266 + t236 * t251 - t240 * t263) * t257) * t278) * t230;
	t1 = [0, t195, 0, 0, t196, 0; 0, (-t200 * t334 - t209 * t324) * t314 + ((t247 * t278 + t254 * t316) * t209 + (-t332 + t338) * t200 + (-t324 * t198 + (-t266 * t316 - t195 * t240 - t204 * t218 - t263 * t278 + (-t204 * t259 - t251 * t278) * t202) * t327 + (-t251 * t316 - t195 * t259 - t204 * t236 - t245 * t278 + (t204 * t240 + t266 * t278) * t202) * t328) * t210) * t206, 0, 0, (-t199 * t334 - t209 * t244) * t314 + (t199 * t338 + t221 * t209 + (-t244 * t198 - t199 * t220 + (-t196 * t240 - t203 * t218 + t237 + (-t203 * t259 - t242) * t202) * t327 + (-t196 * t259 - t203 * t236 - t219 + (t203 * t240 - t260) * t202) * t328) * t210) * t206, 0; 0, 0.2e1 * (-t223 * t234 + t235 * t330) * t335 + (0.2e1 * t235 * t311 + t308 * t223 * t289 + t297 * t331 + (t226 * t287 * t308 - t235 * t213 - t234 * t214 - t297 * t329) * t224) * t215, 0, 0, -t303 * t300 * t313 + (t303 * t220 - ((-qJD(6) * t223 - 0.2e1 * t311) * t289 + (t213 * t289 + (t214 - t315) * t287) * t224) * t300) * t215, t313 + 0.2e1 * (t213 * t224 * t215 + (-t215 * t333 - t224 * t335) * t226) * t226;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end