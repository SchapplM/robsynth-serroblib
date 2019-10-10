% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP3
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
%   Wie in S6RPRRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:55
	% DurationCPUTime: 1.01s
	% Computational Cost: add. (1976->91), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->94)
	t123 = qJ(1) + pkin(10);
	t121 = sin(t123);
	t119 = t121 ^ 2;
	t130 = sin(qJ(3));
	t125 = t130 ^ 2;
	t132 = cos(qJ(3));
	t127 = 0.1e1 / t132 ^ 2;
	t175 = t125 * t127;
	t115 = t119 * t175 + 0.1e1;
	t124 = t130 * t125;
	t126 = 0.1e1 / t132;
	t140 = qJD(3) * (t124 * t127 + t130) * t126;
	t122 = cos(t123);
	t171 = qJD(1) * t122;
	t180 = t121 * t125;
	t145 = t171 * t180;
	t187 = 0.1e1 / t115 ^ 2 * (t119 * t140 + t127 * t145);
	t198 = -0.2e1 * t187;
	t129 = sin(qJ(4));
	t131 = cos(qJ(4));
	t173 = t131 * t132;
	t109 = t121 * t129 + t122 * t173;
	t103 = 0.1e1 / t109;
	t104 = 0.1e1 / t109 ^ 2;
	t174 = t129 * t132;
	t108 = -t121 * t131 + t122 * t174;
	t183 = t104 * t108;
	t142 = -t103 * t129 + t131 * t183;
	t102 = t108 ^ 2;
	t101 = t102 * t104 + 0.1e1;
	t98 = 0.1e1 / t101;
	t197 = t142 * t98;
	t169 = qJD(3) * t121;
	t113 = 0.1e1 / t115;
	t156 = t127 * t169;
	t170 = qJD(1) * t130;
	t157 = t122 * t170;
	t167 = qJD(3) * t132;
	t87 = (-(-t121 * t167 - t157) * t126 + t125 * t156) * t113;
	t149 = t87 - t169;
	t152 = 0.1e1 + t175;
	t196 = t121 * t152;
	t150 = -t121 * t87 + qJD(3);
	t148 = qJD(4) * t132 - qJD(1);
	t168 = qJD(3) * t130;
	t195 = t148 * t129 + t131 * t168;
	t178 = t121 * t130;
	t112 = atan2(-t178, -t132);
	t111 = cos(t112);
	t110 = sin(t112);
	t161 = t110 * t178;
	t96 = -t111 * t132 - t161;
	t93 = 0.1e1 / t96;
	t94 = 0.1e1 / t96 ^ 2;
	t194 = t113 - 0.1e1;
	t120 = t122 ^ 2;
	t162 = t94 * t167;
	t181 = t111 * t130;
	t82 = t150 * t181 + (t149 * t132 - t157) * t110;
	t192 = t82 * t93 * t94;
	t92 = t120 * t125 * t94 + 0.1e1;
	t193 = (-t94 * t145 + (-t125 * t192 + t130 * t162) * t120) / t92 ^ 2;
	t147 = -qJD(1) * t132 + qJD(4);
	t143 = t147 * t131;
	t89 = t121 * t143 - t195 * t122;
	t188 = t103 * t104 * t89;
	t141 = t121 * t174 + t122 * t131;
	t155 = t129 * t168;
	t88 = t141 * qJD(1) - t109 * qJD(4) + t122 * t155;
	t191 = (-t102 * t188 - t88 * t183) / t101 ^ 2;
	t90 = 0.1e1 / t92;
	t190 = t90 * t94;
	t189 = t93 * t90;
	t185 = t122 * t94;
	t182 = t110 * t132;
	t177 = t122 * t129;
	t176 = t122 * t130;
	t172 = qJD(1) * t121;
	t166 = 0.2e1 * t192;
	t165 = -0.2e1 * t191;
	t164 = t93 * t193;
	t163 = t108 * t188;
	t160 = t113 * t125 * t126;
	t158 = t121 * t170;
	t153 = 0.2e1 * t94 * t193;
	t151 = t126 * t198;
	t146 = t121 * t160;
	t144 = t152 * t122;
	t107 = -t121 * t173 + t177;
	t100 = t113 * t196;
	t86 = (t194 * t130 * t110 - t111 * t146) * t122;
	t85 = -t121 * t182 + t181 + (-t111 * t178 + t182) * t100;
	t83 = t196 * t198 + (qJD(1) * t144 + 0.2e1 * t121 * t140) * t113;
	t1 = [t151 * t176 + (qJD(3) * t144 - t126 * t158) * t113, 0, t83, 0, 0, 0; (-t167 * t189 + (0.2e1 * t164 + (qJD(1) * t86 + t82) * t190) * t130) * t121 + (t86 * t153 * t130 + (-t86 * t162 + (t86 * t166 + ((0.2e1 * t130 * t187 - t87 * t146 - t194 * t167) * t110 + (t151 * t180 + t130 * t87 + (t124 * t156 - (t87 - 0.2e1 * t169) * t130) * t113) * t111) * t185) * t130 + (-t93 + (-(t119 - t120) * t111 * t160 + t194 * t161) * t94) * t170) * t90) * t122, 0, (-t172 * t189 + (-0.2e1 * t164 + (-qJD(3) * t85 - t82) * t190) * t122) * t132 + (t85 * t122 * t153 + (-t122 * qJD(3) * t93 - ((-t100 * t171 - t121 * t83) * t111 + (-t150 * t100 - t149) * t110) * t94 * t176 + (t122 * t166 + t94 * t172) * t85 - ((t83 - t171) * t110 + (t149 * t100 + t150) * t111) * t132 * t185) * t90) * t130, 0, 0, 0; 0.2e1 * (t103 * t141 + t107 * t183) * t191 + (0.2e1 * t107 * t163 + (t141 * t89 + t107 * t88 + (-t195 * t121 - t122 * t143) * t108) * t104 + (t147 * t177 + (-t148 * t131 + t155) * t121) * t103) * t98, 0, -t158 * t197 + (t167 * t197 + (t142 * t165 + ((-qJD(4) * t103 - 0.2e1 * t163) * t131 + (-t131 * t88 + (-qJD(4) * t108 + t89) * t129) * t104) * t98) * t130) * t122, t165 + 0.2e1 * (-t104 * t88 * t98 + (-t104 * t191 - t98 * t188) * t108) * t108, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:55
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (2635->97), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
	t158 = sin(qJ(3));
	t153 = t158 ^ 2;
	t159 = cos(qJ(3));
	t155 = 0.1e1 / t159 ^ 2;
	t198 = t153 * t155;
	t151 = qJ(1) + pkin(10);
	t146 = sin(t151);
	t222 = 0.2e1 * t146;
	t221 = t158 * t198;
	t147 = cos(t151);
	t157 = qJ(4) + qJ(5);
	t148 = sin(t157);
	t149 = cos(t157);
	t200 = t149 * t159;
	t134 = t146 * t148 + t147 * t200;
	t150 = qJD(4) + qJD(5);
	t176 = t150 * t159 - qJD(1);
	t193 = qJD(3) * t158;
	t220 = t176 * t148 + t149 * t193;
	t203 = t146 * t158;
	t138 = atan2(-t203, -t159);
	t137 = cos(t138);
	t136 = sin(t138);
	t186 = t136 * t203;
	t124 = -t137 * t159 - t186;
	t121 = 0.1e1 / t124;
	t128 = 0.1e1 / t134;
	t154 = 0.1e1 / t159;
	t122 = 0.1e1 / t124 ^ 2;
	t129 = 0.1e1 / t134 ^ 2;
	t219 = -0.2e1 * t158;
	t144 = t146 ^ 2;
	t142 = t144 * t198 + 0.1e1;
	t140 = 0.1e1 / t142;
	t218 = t140 - 0.1e1;
	t201 = t148 * t159;
	t169 = t146 * t201 + t147 * t149;
	t182 = t148 * t193;
	t112 = t169 * qJD(1) - t134 * t150 + t147 * t182;
	t133 = -t146 * t149 + t147 * t201;
	t127 = t133 ^ 2;
	t120 = t127 * t129 + 0.1e1;
	t208 = t129 * t133;
	t175 = -qJD(1) * t159 + t150;
	t171 = t175 * t149;
	t113 = t146 * t171 - t147 * t220;
	t215 = t113 * t128 * t129;
	t217 = (-t112 * t208 - t127 * t215) / t120 ^ 2;
	t195 = qJD(1) * t158;
	t183 = t147 * t195;
	t192 = qJD(3) * t159;
	t194 = qJD(3) * t146;
	t114 = (-(-t146 * t192 - t183) * t154 + t194 * t198) * t140;
	t206 = t137 * t158;
	t108 = (-t114 * t146 + qJD(3)) * t206 + (-t183 + (t114 - t194) * t159) * t136;
	t216 = t108 * t121 * t122;
	t214 = t114 * t136;
	t213 = t114 * t158;
	t212 = t122 * t147;
	t211 = t122 * t158;
	t168 = qJD(3) * (t158 + t221) * t154;
	t196 = qJD(1) * t147;
	t173 = t146 * t153 * t196;
	t210 = (t144 * t168 + t155 * t173) / t142 ^ 2;
	t180 = 0.1e1 + t198;
	t126 = t180 * t146 * t140;
	t209 = t126 * t146;
	t207 = t136 * t159;
	t145 = t147 ^ 2;
	t205 = t145 * t153;
	t202 = t147 * t148;
	t199 = t153 * t154;
	t197 = qJD(1) * t146;
	t117 = t122 * t205 + 0.1e1;
	t191 = 0.2e1 * (-t205 * t216 + (t145 * t158 * t192 - t173) * t122) / t117 ^ 2;
	t190 = 0.2e1 * t217;
	t189 = 0.2e1 * t216;
	t188 = t133 * t215;
	t187 = t147 * t211;
	t185 = t140 * t199;
	t179 = t158 * t191;
	t178 = t210 * t222;
	t177 = t210 * t219;
	t174 = t146 * t185;
	t172 = t180 * t147;
	t170 = -t128 * t148 + t149 * t208;
	t167 = t170 * t158;
	t132 = -t146 * t200 + t202;
	t118 = 0.1e1 / t120;
	t115 = 0.1e1 / t117;
	t111 = (t218 * t158 * t136 - t137 * t174) * t147;
	t110 = -t146 * t207 + t206 + (-t137 * t203 + t207) * t126;
	t109 = -t180 * t178 + (qJD(1) * t172 + t168 * t222) * t140;
	t105 = -0.2e1 * t217 + 0.2e1 * (-t112 * t118 * t129 + (-t118 * t215 - t129 * t217) * t133) * t133;
	t1 = [t147 * t154 * t177 + (-t146 * t154 * t195 + qJD(3) * t172) * t140, 0, t109, 0, 0, 0; (t121 * t179 + (-t121 * t192 + (qJD(1) * t111 + t108) * t211) * t115) * t146 + (t122 * t179 * t111 + (-((t114 * t174 + t218 * t192 + t177) * t136 + (t178 * t199 - t213 + (t213 + (t219 - t221) * t194) * t140) * t137) * t187 + (-t122 * t192 + t158 * t189) * t111 + (-t121 + ((-t144 + t145) * t137 * t185 + t218 * t186) * t122) * t195) * t115) * t147, 0, (t110 * t211 - t121 * t159) * t147 * t191 + ((-t121 * t197 + (-qJD(3) * t110 - t108) * t212) * t159 + (-t147 * qJD(3) * t121 - (-t109 * t137 * t146 + t136 * t194 + t209 * t214 - t214 + (-qJD(3) * t136 - t137 * t196) * t126) * t187 + (t122 * t197 + t147 * t189) * t110 - ((t109 - t196) * t136 + ((0.1e1 - t209) * qJD(3) + (t126 - t146) * t114) * t137) * t159 * t212) * t158) * t115, 0, 0, 0; (t128 * t169 + t132 * t208) * t190 + (0.2e1 * t132 * t188 + (t132 * t112 + t169 * t113 + (-t146 * t220 - t147 * t171) * t133) * t129 + (t175 * t202 + (-t176 * t149 + t182) * t146) * t128) * t118, 0, -t147 * t167 * t190 + (-t167 * t197 + (t170 * t192 + ((-t128 * t150 - 0.2e1 * t188) * t149 + (-t112 * t149 + (-t133 * t150 + t113) * t148) * t129) * t158) * t147) * t118, t105, t105, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:48:54
	% EndTime: 2019-10-10 01:48:56
	% DurationCPUTime: 1.72s
	% Computational Cost: add. (9974->127), mult. (8378->271), div. (1558->15), fcn. (10537->9), ass. (0->118)
	t192 = qJ(4) + qJ(5);
	t186 = cos(t192);
	t185 = sin(t192);
	t241 = qJ(1) + pkin(10);
	t223 = sin(t241);
	t216 = t223 * t185;
	t184 = cos(t241);
	t194 = cos(qJ(3));
	t250 = t184 * t194;
	t169 = t186 * t250 + t216;
	t163 = 0.1e1 / t169 ^ 2;
	t180 = t184 ^ 2;
	t193 = sin(qJ(3));
	t188 = t193 ^ 2;
	t254 = t180 * t188;
	t232 = t163 * t254;
	t158 = 0.1e1 + t232;
	t212 = qJD(1) * t223;
	t208 = t194 * t212;
	t187 = qJD(4) + qJD(5);
	t214 = t223 * t187;
	t243 = qJD(3) * t193;
	t246 = t187 * t194;
	t148 = (-t208 + t214) * t186 + (-t186 * t243 + (qJD(1) - t246) * t185) * t184;
	t162 = 0.1e1 / t169;
	t261 = t148 * t162 * t163;
	t218 = t254 * t261;
	t242 = qJD(3) * t194;
	t224 = t193 * t242;
	t270 = (-t218 + (-t184 * t188 * t212 + t180 * t224) * t163) / t158 ^ 2;
	t182 = 0.1e1 / t185 ^ 2;
	t248 = t186 * t187;
	t269 = t182 * t248;
	t251 = t184 * t193;
	t211 = t194 * t216;
	t165 = t184 * t186 + t211;
	t226 = t184 * t243;
	t227 = t186 * t246;
	t147 = t165 * qJD(1) - t184 * t227 + (-t214 + t226) * t185;
	t215 = t223 * t186;
	t168 = t185 * t250 - t215;
	t181 = 0.1e1 / t185;
	t189 = 0.1e1 / t193;
	t190 = 0.1e1 / t193 ^ 2;
	t225 = t190 * t242;
	t253 = t181 * t189;
	t268 = (t181 * t225 + t189 * t269) * t168 + t147 * t253;
	t245 = t193 * t185;
	t157 = atan2(-t165, t245);
	t152 = cos(t157);
	t151 = sin(t157);
	t260 = t151 * t165;
	t146 = t152 * t245 - t260;
	t143 = 0.1e1 / t146;
	t144 = 0.1e1 / t146 ^ 2;
	t267 = -0.2e1 * t165;
	t266 = 0.2e1 * t168;
	t160 = t165 ^ 2;
	t252 = t182 * t190;
	t159 = t160 * t252 + 0.1e1;
	t155 = 0.1e1 / t159;
	t247 = t186 * t193;
	t206 = t185 * t242 + t187 * t247;
	t230 = t165 * t252;
	t217 = t193 * t223;
	t209 = qJD(3) * t217;
	t210 = t186 * t212;
	t244 = qJD(1) * t184;
	t249 = t185 * t187;
	t149 = -t185 * t209 - t184 * t249 - t210 + (t185 * t244 + t186 * t214) * t194;
	t233 = t149 * t253;
	t135 = (t206 * t230 - t233) * t155;
	t202 = -t135 * t165 + t206;
	t131 = (-t135 * t245 - t149) * t151 + t202 * t152;
	t145 = t143 * t144;
	t265 = t131 * t145;
	t191 = t189 / t188;
	t228 = t181 * t269;
	t264 = (t149 * t230 + (-t182 * t191 * t242 - t190 * t228) * t160) / t159 ^ 2;
	t263 = t144 * t168;
	t262 = t147 * t144;
	t259 = t151 * t168;
	t258 = t151 * t193;
	t257 = t152 * t165;
	t256 = t152 * t168;
	t255 = t152 * t194;
	t161 = t168 ^ 2;
	t141 = t144 * t161 + 0.1e1;
	t240 = 0.2e1 * (-t161 * t265 - t168 * t262) / t141 ^ 2;
	t239 = -0.2e1 * t264;
	t238 = 0.2e1 * t270;
	t237 = t145 * t266;
	t236 = t189 * t264;
	t235 = t144 * t259;
	t231 = t165 * t253;
	t229 = t181 * t190 * t194;
	t222 = t143 * t240;
	t221 = t144 * t240;
	t220 = t251 * t266;
	t219 = t181 * t236;
	t205 = t165 * t229 + t223;
	t142 = t205 * t155;
	t213 = t223 - t142;
	t167 = -t184 * t185 + t194 * t215;
	t207 = t165 * t182 * t186 - t167 * t181;
	t204 = t163 * t167 * t184 - t223 * t162;
	t153 = 0.1e1 / t158;
	t150 = t169 * qJD(1) - t184 * t248 - t186 * t209 - t187 * t211;
	t139 = 0.1e1 / t141;
	t138 = t207 * t189 * t155;
	t134 = (-t151 + (t152 * t231 + t151) * t155) * t168;
	t133 = -t142 * t257 + (t213 * t258 + t255) * t185;
	t132 = t152 * t247 - t151 * t167 + (-t151 * t245 - t257) * t138;
	t130 = t163 * t220 * t270 + (t220 * t261 + (t147 * t251 + (-t184 * t242 + t193 * t212) * t168) * t163) * t153;
	t129 = t205 * t239 + (t149 * t229 + t244 + (-t227 * t252 + (-0.2e1 * t191 * t194 ^ 2 - t189) * t181 * qJD(3)) * t165) * t155;
	t127 = -0.2e1 * t207 * t236 + (-t207 * t225 + ((-t165 * t187 - t150) * t181 + (t228 * t267 + (t167 * t187 + t149) * t182) * t186) * t189) * t155;
	t126 = (t132 * t263 - t143 * t169) * t240 + (t132 * t262 + t148 * t143 + (t132 * t237 - t144 * t169) * t131 - (t186 * t242 - t187 * t245 - t127 * t165 - t138 * t149 + (-t138 * t245 - t167) * t135) * t144 * t256 - (-t150 + (-t127 * t185 - t135 * t186) * t193 - t202 * t138) * t235) * t139;
	t1 = [t268 * t155 + t219 * t266, 0, t129, t127, t127, 0; t165 * t222 + (-t149 * t143 + (t131 * t165 + t134 * t147) * t144) * t139 + (t134 * t221 + (0.2e1 * t134 * t265 + (t147 * t155 - t147 - (-t135 * t155 * t231 + t239) * t168) * t144 * t151 + (-(t219 * t267 - t135) * t263 + (-(t135 + t233) * t168 + t268 * t165) * t144 * t155) * t152) * t139) * t168, 0, t133 * t168 * t221 + (-(-t129 * t257 + (t135 * t260 - t149 * t152) * t142) * t263 + (-t143 * t251 - (-t142 * t258 + t151 * t217 + t255) * t263) * t248 + (t131 * t237 + t262) * t133) * t139 + (t222 * t251 + ((-t184 * qJD(3) * t143 - (t213 * qJD(3) - t135) * t235) * t194 + (t143 * t212 + (t184 * t131 - (-t129 + t244) * t259 - (t213 * t135 - qJD(3)) * t256) * t144) * t193) * t139) * t185, t126, t126, 0; t204 * t193 * t238 + (-t204 * t242 + ((qJD(1) * t162 + 0.2e1 * t167 * t261) * t184 + (-t223 * t148 - t150 * t184 + t167 * t212) * t163) * t193) * t153, 0, (t162 * t250 + t186 * t232) * t238 + (0.2e1 * t186 * t218 + (t208 + t226) * t162 + ((t148 * t194 + 0.2e1 * t188 * t210) * t184 + (-0.2e1 * t186 * t224 + t188 * t249) * t180) * t163) * t153, t130, t130, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end