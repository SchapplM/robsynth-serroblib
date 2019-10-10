% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR3
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
%   Wie in S6RPRRPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:26:58
	% EndTime: 2019-10-10 01:26:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:26:58
	% EndTime: 2019-10-10 01:26:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:26:58
	% EndTime: 2019-10-10 01:26:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:26:58
	% EndTime: 2019-10-10 01:26:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:26:59
	% EndTime: 2019-10-10 01:27:00
	% DurationCPUTime: 1.02s
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
	% StartTime: 2019-10-10 01:26:59
	% EndTime: 2019-10-10 01:27:00
	% DurationCPUTime: 1.41s
	% Computational Cost: add. (4034->125), mult. (6168->273), div. (1114->15), fcn. (7752->9), ass. (0->111)
	t157 = cos(qJ(4));
	t155 = sin(qJ(4));
	t206 = qJ(1) + pkin(10);
	t185 = sin(t206);
	t180 = t185 * t155;
	t147 = cos(t206);
	t158 = cos(qJ(3));
	t216 = t147 * t158;
	t135 = t157 * t216 + t180;
	t129 = 0.1e1 / t135 ^ 2;
	t146 = t147 ^ 2;
	t156 = sin(qJ(3));
	t151 = t156 ^ 2;
	t218 = t146 * t151;
	t197 = t129 * t218;
	t119 = 0.1e1 + t197;
	t177 = qJD(1) * t185;
	t173 = t158 * t177;
	t176 = t185 * qJD(4);
	t210 = qJD(3) * t156;
	t114 = (-t173 + t176) * t157 + (-t157 * t210 + (-qJD(4) * t158 + qJD(1)) * t155) * t147;
	t128 = 0.1e1 / t135;
	t225 = t114 * t128 * t129;
	t182 = t218 * t225;
	t209 = qJD(3) * t158;
	t231 = (-t182 + (t146 * t156 * t209 - t147 * t151 * t177) * t129) / t119 ^ 2;
	t217 = t147 * t156;
	t131 = t147 * t157 + t158 * t180;
	t172 = t155 * t176;
	t207 = qJD(4) * t157;
	t189 = t147 * t207;
	t193 = t147 * t210;
	t113 = t131 * qJD(1) + t155 * t193 - t158 * t189 - t172;
	t179 = t185 * t157;
	t134 = t155 * t216 - t179;
	t148 = 0.1e1 / t155;
	t149 = 0.1e1 / t155 ^ 2;
	t152 = 0.1e1 / t156;
	t153 = 0.1e1 / t156 ^ 2;
	t192 = t153 * t209;
	t215 = t148 * t152;
	t230 = (t149 * t152 * t207 + t148 * t192) * t134 + t113 * t215;
	t212 = t156 * t155;
	t124 = atan2(-t131, t212);
	t121 = cos(t124);
	t120 = sin(t124);
	t224 = t120 * t131;
	t112 = t121 * t212 - t224;
	t109 = 0.1e1 / t112;
	t110 = 0.1e1 / t112 ^ 2;
	t229 = 0.2e1 * t134;
	t126 = t131 ^ 2;
	t214 = t149 * t153;
	t125 = t126 * t214 + 0.1e1;
	t122 = 0.1e1 / t125;
	t170 = t155 * t209 + t156 * t207;
	t195 = t131 * t214;
	t181 = t156 * t185;
	t174 = qJD(3) * t181;
	t175 = t157 * t177;
	t208 = qJD(4) * t155;
	t211 = qJD(1) * t147;
	t115 = -t155 * t174 - t147 * t208 - t175 + (t155 * t211 + t157 * t176) * t158;
	t198 = t115 * t215;
	t101 = (t170 * t195 - t198) * t122;
	t166 = -t101 * t131 + t170;
	t97 = (-t101 * t212 - t115) * t120 + t166 * t121;
	t228 = t109 * t110 * t97;
	t150 = t148 * t149;
	t154 = t152 / t151;
	t190 = t153 * t207;
	t227 = (t115 * t195 + (-t149 * t154 * t209 - t150 * t190) * t126) / t125 ^ 2;
	t226 = t110 * t134;
	t223 = t120 * t134;
	t222 = t120 * t156;
	t221 = t121 * t131;
	t220 = t121 * t134;
	t219 = t121 * t158;
	t213 = t149 * t157;
	t127 = t134 ^ 2;
	t107 = t110 * t127 + 0.1e1;
	t205 = 0.2e1 / t107 ^ 2 * (-t113 * t226 - t127 * t228);
	t204 = 0.2e1 * t228;
	t203 = 0.2e1 * t231;
	t202 = -0.2e1 * t227;
	t201 = t152 * t227;
	t200 = t110 * t223;
	t196 = t131 * t215;
	t194 = t148 * t153 * t158;
	t191 = t157 * t209;
	t188 = t109 * t205;
	t187 = t110 * t205;
	t186 = t134 * t204;
	t184 = t217 * t229;
	t183 = t148 * t201;
	t169 = t131 * t194 + t185;
	t108 = t169 * t122;
	t178 = t185 - t108;
	t133 = -t147 * t155 + t158 * t179;
	t171 = t131 * t213 - t133 * t148;
	t168 = t129 * t133 * t147 - t185 * t128;
	t117 = 0.1e1 / t119;
	t116 = t135 * qJD(1) - t157 * t174 - t158 * t172 - t189;
	t105 = 0.1e1 / t107;
	t104 = t171 * t152 * t122;
	t100 = (-t120 + (t121 * t196 + t120) * t122) * t134;
	t99 = -t108 * t221 + (t178 * t222 + t219) * t155;
	t98 = t121 * t156 * t157 - t120 * t133 + (-t120 * t212 - t221) * t104;
	t96 = t169 * t202 + (t115 * t194 + t211 + (-t149 * t158 * t190 + (-0.2e1 * t154 * t158 ^ 2 - t152) * t148 * qJD(3)) * t131) * t122;
	t94 = -0.2e1 * t171 * t201 + (-t171 * t192 + (t115 * t213 - t116 * t148 + (t133 * t213 + (-0.2e1 * t150 * t157 ^ 2 - t148) * t131) * qJD(4)) * t152) * t122;
	t1 = [t230 * t122 + t183 * t229, 0, t96, t94, 0, 0; t131 * t188 + (-t115 * t109 + (t100 * t113 + t131 * t97) * t110) * t105 + (t100 * t187 + (t100 * t204 + (t113 * t122 - t113 - (-t101 * t122 * t196 + t202) * t134) * t110 * t120 + (-(-0.2e1 * t131 * t183 - t101) * t226 + (-(t101 + t198) * t134 + t230 * t131) * t110 * t122) * t121) * t105) * t134, 0, t99 * t134 * t187 + (t99 * t186 + (-(-t96 * t221 + (t101 * t224 - t115 * t121) * t108) * t134 + t99 * t113) * t110 + (-t109 * t217 - (-t108 * t222 + t120 * t181 + t219) * t226) * t207) * t105 + (t188 * t217 + ((-t147 * qJD(3) * t109 - (t178 * qJD(3) - t101) * t200) * t158 + (t109 * t177 + (t147 * t97 - (-t96 + t211) * t223 - (t178 * t101 - qJD(3)) * t220) * t110) * t156) * t105) * t155, (-t109 * t135 + t98 * t226) * t205 + (t98 * t186 + t114 * t109 - (-t116 + (-t101 * t157 - t155 * t94) * t156 - t166 * t104) * t200 + (t98 * t113 - t135 * t97 - (t191 - t156 * t208 - t104 * t115 - t131 * t94 + (-t104 * t212 - t133) * t101) * t220) * t110) * t105, 0, 0; t168 * t156 * t203 + (-t168 * t209 + ((qJD(1) * t128 + 0.2e1 * t133 * t225) * t147 + (-t185 * t114 - t116 * t147 + t133 * t177) * t129) * t156) * t117, 0, (t128 * t216 + t157 * t197) * t203 + (0.2e1 * t157 * t182 + (t173 + t193) * t128 + ((t114 * t158 + 0.2e1 * t151 * t175) * t147 + (t151 * t208 - 0.2e1 * t156 * t191) * t146) * t129) * t117, t129 * t184 * t231 + (t184 * t225 + (t113 * t217 + (-t147 * t209 + t156 * t177) * t134) * t129) * t117, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:26:59
	% EndTime: 2019-10-10 01:27:00
	% DurationCPUTime: 1.47s
	% Computational Cost: add. (2901->118), mult. (4276->251), div. (493->12), fcn. (5073->11), ass. (0->114)
	t188 = sin(qJ(3));
	t182 = t188 ^ 2;
	t191 = cos(qJ(3));
	t184 = 0.1e1 / t191 ^ 2;
	t239 = t182 * t184;
	t180 = qJ(1) + pkin(10);
	t178 = sin(t180);
	t176 = t178 ^ 2;
	t172 = t176 * t239 + 0.1e1;
	t183 = 0.1e1 / t191;
	t265 = t188 * t239;
	t199 = qJD(3) * (t188 + t265) * t183;
	t179 = cos(t180);
	t235 = qJD(1) * t179;
	t243 = t178 * t182;
	t208 = t235 * t243;
	t250 = (t176 * t199 + t184 * t208) / t172 ^ 2;
	t267 = -0.2e1 * t250;
	t215 = 0.1e1 + t239;
	t266 = t178 * t215;
	t190 = cos(qJ(4));
	t264 = (-qJD(4) + qJD(6)) * t190;
	t187 = sin(qJ(4));
	t230 = qJD(4) * t187;
	t263 = -qJD(6) * t187 + t230;
	t211 = qJD(4) * t191 - qJD(1);
	t232 = qJD(3) * t188;
	t262 = t211 * t187 + t190 * t232;
	t238 = t187 * t191;
	t200 = t178 * t238 + t179 * t190;
	t217 = t187 * t232;
	t237 = t190 * t191;
	t220 = t179 * t237;
	t137 = t200 * qJD(1) - qJD(4) * t220 - t178 * t230 + t179 * t217;
	t210 = -qJD(1) * t191 + qJD(4);
	t201 = t210 * t190;
	t138 = t178 * t201 - t262 * t179;
	t186 = sin(qJ(6));
	t189 = cos(qJ(6));
	t165 = -t178 * t190 + t179 * t238;
	t166 = t178 * t187 + t220;
	t204 = t165 * t189 - t166 * t186;
	t130 = t204 * qJD(6) - t137 * t186 + t138 * t189;
	t153 = t165 * t186 + t166 * t189;
	t145 = 0.1e1 / t153;
	t202 = t186 * t187 + t189 * t190;
	t203 = t186 * t190 - t187 * t189;
	t146 = 0.1e1 / t153 ^ 2;
	t252 = t146 * t204;
	t261 = t203 * t145 + t202 * t252;
	t242 = t178 * t188;
	t171 = atan2(t242, t191);
	t168 = cos(t171);
	t167 = sin(t171);
	t222 = t167 * t242;
	t159 = t168 * t191 + t222;
	t156 = 0.1e1 / t159;
	t157 = 0.1e1 / t159 ^ 2;
	t260 = 0.2e1 * t188;
	t169 = 0.1e1 / t172;
	t259 = t169 - 0.1e1;
	t129 = t153 * qJD(6) + t137 * t189 + t138 * t186;
	t144 = t204 ^ 2;
	t135 = t144 * t146 + 0.1e1;
	t147 = t145 * t146;
	t255 = t130 * t147;
	t258 = (-t129 * t252 - t144 * t255) / t135 ^ 2;
	t177 = t179 ^ 2;
	t244 = t177 * t182;
	t143 = t157 * t244 + 0.1e1;
	t231 = qJD(3) * t191;
	t234 = qJD(1) * t188;
	t218 = t179 * t234;
	t233 = qJD(3) * t178;
	t136 = ((t178 * t231 + t218) * t183 + t233 * t239) * t169;
	t245 = t168 * t188;
	t127 = (t136 * t178 - qJD(3)) * t245 + (t218 + (-t136 + t233) * t191) * t167;
	t256 = t127 * t156 * t157;
	t257 = (-t244 * t256 + (t177 * t188 * t231 - t208) * t157) / t143 ^ 2;
	t254 = t136 * t167;
	t253 = t136 * t188;
	t240 = t179 * t188;
	t161 = t202 * t240;
	t251 = t146 * t161;
	t155 = t169 * t266;
	t249 = t155 * t178;
	t248 = t157 * t179;
	t247 = t157 * t188;
	t246 = t167 * t191;
	t241 = t179 * t187;
	t236 = qJD(1) * t178;
	t226 = 0.2e1 * t258;
	t225 = -0.2e1 * t256;
	t224 = -0.2e1 * t147 * t204;
	t223 = t157 * t240;
	t221 = t169 * t182 * t183;
	t219 = t178 * t234;
	t214 = -0.2e1 * t188 * t257;
	t213 = t130 * t224;
	t212 = t183 * t267;
	t209 = t178 * t221;
	t207 = t215 * t179;
	t164 = -t178 * t237 + t241;
	t205 = -t164 * t186 - t189 * t200;
	t149 = t164 * t189 - t186 * t200;
	t160 = t203 * t240;
	t141 = 0.1e1 / t143;
	t140 = t262 * t178 + t179 * t201;
	t139 = t210 * t241 + (-t211 * t190 + t217) * t178;
	t133 = 0.1e1 / t135;
	t132 = (-t259 * t188 * t167 + t168 * t209) * t179;
	t131 = t178 * t246 - t245 + (t168 * t242 - t246) * t155;
	t128 = t266 * t267 + (qJD(1) * t207 + 0.2e1 * t178 * t199) * t169;
	t1 = [t212 * t240 + (qJD(3) * t207 - t183 * t219) * t169, 0, t128, 0, 0, 0; (t156 * t214 + (t156 * t231 + (-qJD(1) * t132 - t127) * t247) * t141) * t178 + (t157 * t214 * t132 + (((-t136 * t209 - t259 * t231 + t250 * t260) * t167 + (t212 * t243 + t253 + (-t253 + (t260 + t265) * t233) * t169) * t168) * t223 + (t157 * t231 + t188 * t225) * t132 + (t156 + ((-t176 + t177) * t168 * t221 + t259 * t222) * t157) * t234) * t141) * t179, 0, 0.2e1 * (-t131 * t247 + t156 * t191) * t179 * t257 + ((t156 * t236 + (qJD(3) * t131 + t127) * t248) * t191 + (t179 * qJD(3) * t156 + (t128 * t168 * t178 - t167 * t233 - t249 * t254 + t254 + (qJD(3) * t167 + t168 * t235) * t155) * t223 + (-t157 * t236 + t179 * t225) * t131 + ((-t128 + t235) * t167 + ((-0.1e1 + t249) * qJD(3) + (-t155 + t178) * t136) * t168) * t191 * t248) * t188) * t141, 0, 0, 0; (t145 * t205 - t149 * t252) * t226 + ((t149 * qJD(6) - t139 * t189 + t140 * t186) * t145 + t149 * t213 + (t205 * t130 + (t205 * qJD(6) + t139 * t186 + t140 * t189) * t204 - t149 * t129) * t146) * t133, 0, (t145 * t160 + t204 * t251) * t226 + (t129 * t251 + (t146 * t160 - t161 * t224) * t130 + t261 * t219 + (-t261 * t231 + ((-t264 * t145 + t263 * t252) * t189 + (t263 * t145 + t264 * t252) * t186) * t188) * t179) * t133, (t145 * t153 + t204 * t252) * t226 + (-t130 * t145 - t204 * t213 + (0.2e1 * t204 * t129 + t153 * t130) * t146) * t133, 0, -0.2e1 * t258 - 0.2e1 * (t129 * t133 * t146 - (-t133 * t255 - t146 * t258) * t204) * t204;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end