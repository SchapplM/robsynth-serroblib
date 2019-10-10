% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP4
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
%   Wie in S6RRRPRP4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:40
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:27
	% DurationCPUTime: 0.92s
	% Computational Cost: add. (3078->72), mult. (2858->158), div. (686->14), fcn. (3330->7), ass. (0->75)
	t122 = sin(qJ(1));
	t115 = t122 ^ 2;
	t121 = qJ(2) + qJ(3);
	t112 = sin(t121);
	t107 = t112 ^ 2;
	t113 = cos(t121);
	t110 = 0.1e1 / t113 ^ 2;
	t156 = t107 * t110;
	t102 = t115 * t156 + 0.1e1;
	t106 = t112 * t107;
	t108 = t113 ^ 2;
	t109 = 0.1e1 / t113;
	t114 = qJD(2) + qJD(3);
	t155 = t109 * t112;
	t131 = t114 * (t106 * t109 / t108 + t155);
	t123 = cos(qJ(1));
	t147 = qJD(1) * t123;
	t139 = t122 * t147;
	t164 = 0.1e1 / t102 ^ 2 * (t115 * t131 + t139 * t156);
	t172 = -0.2e1 * t164;
	t100 = 0.1e1 / t102;
	t137 = 0.1e1 + t156;
	t170 = t122 * t137;
	t95 = t100 * t170;
	t171 = t122 * t95 - 0.1e1;
	t151 = 0.1e1 / t122 * t123;
	t120 = t123 ^ 2;
	t169 = qJD(1) * (0.1e1 / t115 * t120 + 0.1e1) * t151;
	t149 = t122 * t112;
	t99 = atan2(-t149, -t113);
	t97 = sin(t99);
	t143 = t97 * t149;
	t98 = cos(t99);
	t94 = -t113 * t98 - t143;
	t91 = 0.1e1 / t94;
	t92 = 0.1e1 / t94 ^ 2;
	t153 = t113 * t114;
	t135 = t112 * t120 * t153;
	t152 = t114 * t122;
	t162 = t113 * t97;
	t140 = t110 * t152;
	t86 = (-(-t112 * t147 - t113 * t152) * t109 + t107 * t140) * t100;
	t81 = (t86 - t152) * t162 + (-t97 * t147 + (-t122 * t86 + t114) * t98) * t112;
	t167 = t81 * t91 * t92;
	t89 = t107 * t120 * t92 + 0.1e1;
	t168 = (t92 * t135 + (-t120 * t167 - t92 * t139) * t107) / t89 ^ 2;
	t87 = 0.1e1 / t89;
	t165 = t87 * t92;
	t163 = t112 * t97;
	t161 = t114 * t95;
	t159 = t123 * t92;
	t158 = t98 * t112;
	t157 = t107 * t109;
	t154 = t112 * t123;
	t117 = 0.1e1 / t122 ^ 2;
	t150 = t117 * t120;
	t148 = qJD(1) * t122;
	t146 = 0.2e1 * t167;
	t105 = t108 * t150 + 0.1e1;
	t145 = 0.2e1 / t105 ^ 2 * (-t108 * t169 - t117 * t135);
	t144 = t91 * t168;
	t142 = t87 * t153;
	t141 = t122 * t157;
	t138 = 0.2e1 * t92 * t168;
	t136 = 0.1e1 + t150;
	t134 = t137 * t123;
	t133 = t136 * t112;
	t130 = -t98 * t141 + t163;
	t103 = 0.1e1 / t105;
	t85 = (t130 * t100 - t163) * t123;
	t84 = t112 * t145 * t151 + (qJD(1) * t133 - t151 * t153) * t103;
	t83 = (-t122 + t95) * t162 - t171 * t158;
	t82 = t170 * t172 + (qJD(1) * t134 + 0.2e1 * t122 * t131) * t100;
	t79 = (-t91 * t87 * t148 + (-0.2e1 * t144 + (-t114 * t83 - t81) * t165) * t123) * t113 + (t83 * t123 * t138 + (-t123 * t114 * t91 - ((-t122 * t82 - t147 * t95) * t98 + (t171 * t86 + t152 - t161) * t97) * t92 * t154 + (t123 * t146 + t92 * t148) * t83 - ((t82 - t147) * t97 + (t86 * t95 + t114 + (-t86 - t161) * t122) * t98) * t113 * t159) * t87) * t112;
	t1 = [t109 * t154 * t172 + (t114 * t134 - t148 * t155) * t100, t82, t82, 0, 0, 0; (-t91 * t142 + (0.2e1 * t144 + (qJD(1) * t85 + t81) * t165) * t112) * t122 + (-t85 * t92 * t142 + (t85 * t138 + (t85 * t146 + (t86 * t158 + t97 * t153 + 0.2e1 * t130 * t164 + ((-t86 * t141 - t153) * t97 + (t106 * t140 - (t86 - 0.2e1 * t152) * t112) * t98) * t100) * t159) * t87) * t112 + (-t91 + (-t143 + (t143 - (t115 - t120) * t98 * t157) * t100) * t92) * t112 * t87 * qJD(1)) * t123, t79, t79, 0, 0, 0; t136 * t113 * t145 + (0.2e1 * t113 * t169 + t114 * t133) * t103, t84, t84, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:27
	% DurationCPUTime: 1.07s
	% Computational Cost: add. (3360->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->95)
	t173 = sin(qJ(1));
	t171 = qJ(2) + qJ(3);
	t166 = sin(t171);
	t162 = 0.1e1 / t166 ^ 2;
	t167 = cos(t171);
	t165 = t167 ^ 2;
	t219 = t162 * t165;
	t194 = 0.1e1 + t219;
	t234 = t173 * t194;
	t169 = t173 ^ 2;
	t159 = t169 * t219 + 0.1e1;
	t157 = 0.1e1 / t159;
	t161 = 0.1e1 / t166;
	t175 = cos(qJ(1));
	t206 = qJD(1) * t175;
	t195 = t167 * t206;
	t168 = qJD(2) + qJD(3);
	t214 = t168 * t173;
	t197 = t162 * t214;
	t131 = ((t166 * t214 - t195) * t161 + t165 * t197) * t157;
	t233 = -t131 + t214;
	t190 = qJD(1) * t166 + qJD(5);
	t213 = t168 * t175;
	t232 = -t167 * t213 + t190 * t173;
	t212 = t173 * t167;
	t156 = atan2(-t212, t166);
	t155 = cos(t156);
	t154 = sin(t156);
	t199 = t154 * t212;
	t141 = t155 * t166 - t199;
	t138 = 0.1e1 / t141;
	t172 = sin(qJ(5));
	t209 = t175 * t172;
	t174 = cos(qJ(5));
	t210 = t173 * t174;
	t151 = t166 * t209 + t210;
	t147 = 0.1e1 / t151;
	t139 = 0.1e1 / t141 ^ 2;
	t148 = 0.1e1 / t151 ^ 2;
	t231 = t157 - 0.1e1;
	t221 = t155 * t167;
	t126 = (-t131 * t173 + t168) * t221 + (t233 * t166 - t195) * t154;
	t230 = t126 * t138 * t139;
	t191 = qJD(5) * t166 + qJD(1);
	t186 = t191 * t175;
	t135 = t172 * t186 + t232 * t174;
	t208 = t175 * t174;
	t211 = t173 * t172;
	t150 = -t166 * t208 + t211;
	t146 = t150 ^ 2;
	t145 = t146 * t148 + 0.1e1;
	t224 = t148 * t150;
	t136 = -t232 * t172 + t174 * t186;
	t228 = t136 * t147 * t148;
	t229 = (t135 * t224 - t146 * t228) / t145 ^ 2;
	t164 = t167 * t165;
	t220 = t161 * t167;
	t184 = t168 * (-t161 * t162 * t164 - t220);
	t217 = t165 * t173;
	t188 = t206 * t217;
	t227 = (t162 * t188 + t169 * t184) / t159 ^ 2;
	t226 = t139 * t167;
	t225 = t139 * t175;
	t223 = t150 * t172;
	t222 = t154 * t173;
	t170 = t175 ^ 2;
	t218 = t165 * t170;
	t216 = t166 * t168;
	t215 = t167 * t168;
	t207 = qJD(1) * t173;
	t134 = t139 * t218 + 0.1e1;
	t205 = 0.2e1 * (-t218 * t230 + (-t166 * t170 * t215 - t188) * t139) / t134 ^ 2;
	t204 = 0.2e1 * t230;
	t203 = 0.2e1 * t229;
	t202 = -0.2e1 * t227;
	t201 = t167 * t227;
	t200 = t167 * t225;
	t198 = t161 * t217;
	t193 = t167 * t205;
	t192 = 0.2e1 * t150 * t228;
	t189 = t155 * t157 * t161 * t165;
	t187 = t194 * t175;
	t185 = t147 * t174 + t148 * t223;
	t183 = t185 * t175;
	t153 = -t166 * t211 + t208;
	t152 = t166 * t210 + t209;
	t143 = 0.1e1 / t145;
	t142 = t157 * t234;
	t132 = 0.1e1 / t134;
	t130 = (t231 * t167 * t154 + t173 * t189) * t175;
	t128 = t166 * t222 + t221 + (-t154 * t166 - t155 * t212) * t142;
	t127 = t202 * t234 + (qJD(1) * t187 + 0.2e1 * t173 * t184) * t157;
	t124 = t167 * t183 * t203 + (t183 * t216 + (t185 * t207 + ((qJD(5) * t147 + t192) * t172 + (-t135 * t172 + (-qJD(5) * t150 + t136) * t174) * t148) * t175) * t167) * t143;
	t123 = (t128 * t226 + t138 * t166) * t175 * t205 + ((t138 * t207 + (t128 * t168 + t126) * t225) * t166 + (-t138 * t213 - (-t127 * t155 * t173 + t233 * t154 + (t131 * t222 - t154 * t168 - t155 * t206) * t142) * t200 + (t139 * t207 + t175 * t204) * t128 - ((-t127 + t206) * t154 + ((t142 * t173 - 0.1e1) * t168 + (-t142 + t173) * t131) * t155) * t166 * t225) * t167) * t132;
	t1 = [0.2e1 * t175 * t161 * t201 + (t168 * t187 + t207 * t220) * t157, t127, t127, 0, 0, 0; (t138 * t193 + (t138 * t216 + (qJD(1) * t130 + t126) * t226) * t132) * t173 + (t139 * t193 * t130 + (-((-0.2e1 * t201 + t216 + (-t131 * t198 - t216) * t157) * t154 + (t198 * t202 - t131 * t167 + (-t164 * t197 + (t131 - 0.2e1 * t214) * t167) * t157) * t155) * t200 + (t139 * t216 + t167 * t204) * t130 + (-t138 + ((t169 - t170) * t189 + t231 * t199) * t139) * t167 * qJD(1)) * t132) * t175, t123, t123, 0, 0, 0; (-t147 * t152 + t153 * t224) * t203 + (t153 * t192 + (-t153 * t135 - t152 * t136 + t191 * t150 * t210 - (-t168 * t212 - t190 * t175) * t223) * t148 + (t190 * t208 + (-t191 * t172 + t174 * t215) * t173) * t147) * t143, t124, t124, 0, -0.2e1 * t229 + 0.2e1 * (t135 * t148 * t143 + (-t143 * t228 - t148 * t229) * t150) * t150, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:40:26
	% EndTime: 2019-10-10 11:40:28
	% DurationCPUTime: 1.54s
	% Computational Cost: add. (5529->122), mult. (8382->269), div. (1515->15), fcn. (10508->9), ass. (0->120)
	t194 = qJ(2) + qJ(3);
	t187 = sin(t194);
	t195 = sin(qJ(5));
	t198 = cos(qJ(1));
	t243 = t198 * t195;
	t196 = sin(qJ(1));
	t197 = cos(qJ(5));
	t245 = t196 * t197;
	t174 = t187 * t243 + t245;
	t170 = 0.1e1 / t174 ^ 2;
	t188 = cos(t194);
	t183 = t188 ^ 2;
	t193 = t198 ^ 2;
	t256 = t183 * t193;
	t222 = t170 * t256;
	t166 = 0.1e1 + t222;
	t216 = qJD(5) * t187 + qJD(1);
	t212 = t216 * t197;
	t215 = qJD(1) * t187 + qJD(5);
	t189 = qJD(2) + qJD(3);
	t248 = t189 * t198;
	t223 = t188 * t248;
	t158 = t198 * t212 + (-t215 * t196 + t223) * t195;
	t169 = 0.1e1 / t174;
	t263 = t158 * t169 * t170;
	t214 = t256 * t263;
	t239 = qJD(1) * t198;
	t220 = t196 * t239;
	t225 = t188 * t189 * t193;
	t272 = (-t214 + (-t183 * t220 - t187 * t225) * t170) / t166 ^ 2;
	t184 = 0.1e1 / t188;
	t175 = t187 * t245 + t243;
	t242 = t198 * t197;
	t246 = t196 * t195;
	t176 = -t187 * t246 + t242;
	t190 = 0.1e1 / t197;
	t191 = 0.1e1 / t197 ^ 2;
	t247 = t191 * t195;
	t209 = t175 * t247 + t176 * t190;
	t271 = t184 * t209;
	t250 = t188 * t198;
	t251 = t188 * t197;
	t165 = atan2(t175, t251);
	t160 = cos(t165);
	t159 = sin(t165);
	t261 = t159 * t175;
	t154 = t160 * t251 + t261;
	t151 = 0.1e1 / t154;
	t152 = 0.1e1 / t154 ^ 2;
	t221 = t187 * t242;
	t173 = -t221 + t246;
	t270 = 0.2e1 * t173;
	t269 = -0.2e1 * t195;
	t168 = t173 ^ 2;
	t149 = t168 * t152 + 0.1e1;
	t157 = t175 * qJD(1) + t174 * qJD(5) - t197 * t223;
	t264 = t157 * t152;
	t172 = t175 ^ 2;
	t185 = 0.1e1 / t188 ^ 2;
	t254 = t185 * t191;
	t167 = t172 * t254 + 0.1e1;
	t163 = 0.1e1 / t167;
	t238 = qJD(5) * t195;
	t253 = t187 * t189;
	t208 = -t188 * t238 - t197 * t253;
	t228 = t175 * t254;
	t252 = t188 * t196;
	t224 = t189 * t252;
	t237 = qJD(5) * t197;
	t155 = -qJD(1) * t221 - t197 * t224 - t198 * t237 + t216 * t246;
	t255 = t184 * t190;
	t230 = t155 * t255;
	t143 = (-t208 * t228 - t230) * t163;
	t207 = -t143 * t175 - t208;
	t139 = (-t143 * t251 - t155) * t159 - t207 * t160;
	t153 = t151 * t152;
	t267 = t139 * t153;
	t268 = (-t168 * t267 + t173 * t264) / t149 ^ 2;
	t186 = t184 / t183;
	t192 = t190 * t191;
	t266 = (-t155 * t228 + (t185 * t192 * t238 + t186 * t191 * t253) * t172) / t167 ^ 2;
	t265 = t152 * t173;
	t262 = t159 * t173;
	t260 = t159 * t188;
	t259 = t160 * t173;
	t258 = t160 * t175;
	t257 = t160 * t187;
	t249 = t189 * t190;
	t244 = t198 * t151;
	t227 = t185 * t187 * t190;
	t211 = t175 * t227 + t196;
	t150 = t211 * t163;
	t241 = t150 - t196;
	t240 = qJD(1) * t196;
	t236 = 0.2e1 * t268;
	t235 = 0.2e1 * t272;
	t234 = -0.2e1 * t266;
	t233 = t153 * t270;
	t232 = t151 * t268;
	t231 = t152 * t262;
	t229 = t175 * t255;
	t226 = t187 * t249;
	t219 = t191 * t238;
	t218 = t152 * t236;
	t217 = t250 * t270;
	t213 = 0.2e1 * t255 * t266;
	t210 = -t170 * t176 * t198 - t169 * t196;
	t206 = t157 * t255 - (-t184 * t219 - t185 * t226) * t173;
	t161 = 0.1e1 / t166;
	t156 = t196 * t212 + (t215 * t198 + t224) * t195;
	t147 = 0.1e1 / t149;
	t146 = t163 * t271;
	t142 = (-t159 + (-t160 * t229 + t159) * t163) * t173;
	t141 = t150 * t258 + (-t241 * t260 - t257) * t197;
	t140 = -t160 * t188 * t195 + t159 * t176 + (-t159 * t251 + t258) * t146;
	t138 = t211 * t234 + (-t155 * t227 + t239 + (t184 * t249 + (t185 * t219 + 0.2e1 * t186 * t226) * t187) * t175) * t163;
	t136 = t234 * t271 + (t209 * t185 * t253 + (-t155 * t247 - t156 * t190 + (t176 * t247 + (0.2e1 * t192 * t195 ^ 2 + t190) * t175) * qJD(5)) * t184) * t163;
	t135 = (-t169 * t187 * t198 - t195 * t222) * t235 + (t214 * t269 + (-t187 * t240 + t223) * t169 + ((-t158 * t198 + t225 * t269) * t187 + (t193 * t237 + t220 * t269) * t183) * t170) * t161;
	t134 = t141 * t173 * t218 + (-(t138 * t258 + (-t143 * t261 - t155 * t160) * t150) * t265 + (t139 * t233 - t264) * t141 + (t188 * t244 - (t150 * t260 - t159 * t252 + t257) * t265) * t238) * t147 + (0.2e1 * t232 * t250 + ((t189 * t244 - (t241 * t189 + t143) * t231) * t187 + (t151 * t240 + (t198 * t139 - (-t138 + t239) * t262 - (-t241 * t143 - t189) * t259) * t152) * t188) * t147) * t197;
	t1 = [-t206 * t163 + t173 * t213, t138, t138, 0, t136, 0; -0.2e1 * t175 * t232 + (-t155 * t151 + (-t139 * t175 - t142 * t157) * t152) * t147 + (t142 * t218 + (0.2e1 * t142 * t267 + (-t157 * t163 + t157 - (t143 * t163 * t229 + t234) * t173) * t152 * t159 + (-(t175 * t213 - t143) * t265 + (-(t143 + t230) * t173 + t206 * t175) * t152 * t163) * t160) * t147) * t173, t134, t134, 0, (t140 * t265 - t151 * t174) * t236 + (-t140 * t264 + t158 * t151 + (t140 * t233 - t174 * t152) * t139 - (-t188 * t237 + t195 * t253 + t136 * t175 - t146 * t155 + (-t146 * t251 + t176) * t143) * t152 * t259 - (-t156 + (-t136 * t197 + t143 * t195) * t188 + t207 * t146) * t231) * t147, 0; t210 * t188 * t235 + (t210 * t253 + ((qJD(1) * t169 - 0.2e1 * t176 * t263) * t198 + (-t156 * t198 + (-qJD(1) * t176 - t158) * t196) * t170) * t188) * t161, t135, t135, 0, t170 * t217 * t272 + (t217 * t263 + (-t157 * t250 + (t187 * t248 + t188 * t240) * t173) * t170) * t161, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end