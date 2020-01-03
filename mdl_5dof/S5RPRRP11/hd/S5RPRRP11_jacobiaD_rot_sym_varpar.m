% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRP11
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
%   Wie in S5RPRRP11_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRP11_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP11_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:55:25
	% EndTime: 2019-12-31 18:55:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:55:25
	% EndTime: 2019-12-31 18:55:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:55:25
	% EndTime: 2019-12-31 18:55:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:55:25
	% EndTime: 2019-12-31 18:55:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:55:25
	% EndTime: 2019-12-31 18:55:26
	% DurationCPUTime: 0.72s
	% Computational Cost: add. (2270->94), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->93)
	t136 = sin(qJ(1));
	t133 = t136 ^ 2;
	t132 = pkin(8) + qJ(3);
	t130 = sin(t132);
	t126 = t130 ^ 2;
	t131 = cos(t132);
	t128 = 0.1e1 / t131 ^ 2;
	t185 = t126 * t128;
	t121 = t133 * t185 + 0.1e1;
	t125 = t130 * t126;
	t127 = 0.1e1 / t131;
	t182 = t127 * t130;
	t146 = qJD(3) * (t125 * t127 * t128 + t182);
	t138 = cos(qJ(1));
	t174 = qJD(1) * t138;
	t183 = t126 * t136;
	t151 = t174 * t183;
	t191 = (t128 * t151 + t133 * t146) / t121 ^ 2;
	t201 = -0.2e1 * t191;
	t158 = 0.1e1 + t185;
	t200 = t136 * t158;
	t137 = cos(qJ(4));
	t176 = t137 * t138;
	t135 = sin(qJ(4));
	t178 = t136 * t135;
	t117 = t131 * t176 + t178;
	t179 = t136 * t130;
	t118 = atan2(-t179, -t131);
	t113 = cos(t118);
	t112 = sin(t118);
	t164 = t112 * t179;
	t102 = -t113 * t131 - t164;
	t99 = 0.1e1 / t102;
	t109 = 0.1e1 / t117;
	t100 = 0.1e1 / t102 ^ 2;
	t110 = 0.1e1 / t117 ^ 2;
	t119 = 0.1e1 / t121;
	t199 = t119 - 0.1e1;
	t134 = t138 ^ 2;
	t173 = qJD(3) * t131;
	t184 = t126 * t134;
	t172 = qJD(3) * t136;
	t160 = t128 * t172;
	t161 = t130 * t174;
	t93 = (-(-t131 * t172 - t161) * t127 + t126 * t160) * t119;
	t155 = t93 - t172;
	t156 = -t136 * t93 + qJD(3);
	t187 = t113 * t130;
	t88 = t156 * t187 + (t155 * t131 - t161) * t112;
	t196 = t99 * t100 * t88;
	t96 = t100 * t184 + 0.1e1;
	t198 = (-t184 * t196 + (t130 * t134 * t173 - t151) * t100) / t96 ^ 2;
	t94 = 0.1e1 / t96;
	t197 = t100 * t94;
	t177 = t136 * t137;
	t180 = t135 * t138;
	t116 = t131 * t180 - t177;
	t108 = t116 ^ 2;
	t107 = t108 * t110 + 0.1e1;
	t189 = t110 * t116;
	t153 = -qJD(1) * t131 + qJD(4);
	t154 = qJD(4) * t131 - qJD(1);
	t171 = qJD(3) * t138;
	t159 = t130 * t171;
	t98 = -t154 * t180 + (t153 * t136 - t159) * t137;
	t194 = t109 * t110 * t98;
	t147 = t131 * t178 + t176;
	t97 = t147 * qJD(1) - t117 * qJD(4) + t135 * t159;
	t195 = 0.1e1 / t107 ^ 2 * (-t108 * t194 - t97 * t189);
	t190 = t109 * t135;
	t188 = t112 * t131;
	t186 = t116 * t137;
	t181 = t130 * t138;
	t175 = qJD(1) * t136;
	t170 = 0.2e1 * t198;
	t169 = 0.2e1 * t196;
	t168 = -0.2e1 * t195;
	t167 = t99 * t198;
	t166 = t116 * t194;
	t165 = t94 * t173;
	t163 = t119 * t126 * t127;
	t157 = t127 * t201;
	t152 = t136 * t163;
	t150 = t158 * t138;
	t149 = t153 * t138;
	t148 = t110 * t186 - t190;
	t115 = -t131 * t177 + t180;
	t105 = 0.1e1 / t107;
	t104 = t119 * t200;
	t92 = (t199 * t130 * t112 - t113 * t152) * t138;
	t90 = -t136 * t188 + t187 + (-t113 * t179 + t188) * t104;
	t89 = t200 * t201 + (qJD(1) * t150 + 0.2e1 * t136 * t146) * t119;
	t1 = [t157 * t181 + (qJD(3) * t150 - t175 * t182) * t119, 0, t89, 0, 0; (-t99 * t165 + (0.2e1 * t167 + (qJD(1) * t92 + t88) * t197) * t130) * t136 + ((-t92 * t165 + (t92 * t170 + ((0.2e1 * t130 * t191 - t93 * t152 - t199 * t173) * t112 + (t157 * t183 + t130 * t93 + (t125 * t160 - (t93 - 0.2e1 * t172) * t130) * t119) * t113) * t94 * t138) * t130) * t100 + (t92 * t169 + (-t99 + ((-t133 + t134) * t113 * t163 + t199 * t164) * t100) * qJD(1)) * t130 * t94) * t138, 0, (-t99 * t94 * t175 + (-0.2e1 * t167 + (-qJD(3) * t90 - t88) * t197) * t138) * t131 + (((-qJD(3) * t99 + t90 * t169) * t138 + (t90 * t175 + (-(-t104 * t174 - t136 * t89) * t113 - ((t104 * t136 - 0.1e1) * t93 + (-t104 + t136) * qJD(3)) * t112) * t181) * t100) * t94 + (t90 * t170 - ((t89 - t174) * t112 + (t155 * t104 + t156) * t113) * t94 * t131) * t100 * t138) * t130, 0, 0; 0.2e1 * (t109 * t147 + t115 * t189) * t195 + (0.2e1 * t115 * t166 - t154 * t109 * t177 + (t130 * t172 + t149) * t190 + (t147 * t98 + t115 * t97 - t149 * t186 - (qJD(3) * t130 * t137 + t154 * t135) * t116 * t136) * t110) * t105, 0, t148 * t168 * t181 + (t148 * t131 * t171 + (-t148 * t175 + ((-qJD(4) * t109 - 0.2e1 * t166) * t137 + (-t137 * t97 + (-qJD(4) * t116 + t98) * t135) * t110) * t138) * t130) * t105, t168 + 0.2e1 * (-t105 * t110 * t97 + (-t105 * t194 - t110 * t195) * t116) * t116, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:55:25
	% EndTime: 2019-12-31 18:55:26
	% DurationCPUTime: 0.96s
	% Computational Cost: add. (3877->124), mult. (6168->272), div. (1114->15), fcn. (7752->9), ass. (0->115)
	t159 = pkin(8) + qJ(3);
	t158 = cos(t159);
	t164 = sin(qJ(4));
	t238 = sin(qJ(1));
	t198 = t238 * t164;
	t165 = cos(qJ(4));
	t166 = cos(qJ(1));
	t218 = t166 * t165;
	t142 = t158 * t218 + t198;
	t136 = 0.1e1 / t142 ^ 2;
	t157 = sin(t159);
	t153 = t157 ^ 2;
	t163 = t166 ^ 2;
	t226 = t153 * t163;
	t203 = t136 * t226;
	t131 = 0.1e1 + t203;
	t190 = qJD(1) * t238;
	t215 = qJD(3) * t166;
	t194 = t157 * t215;
	t176 = t158 * t190 + t194;
	t189 = t238 * qJD(4);
	t219 = t166 * t164;
	t121 = (-qJD(4) * t158 + qJD(1)) * t219 + (t189 - t176) * t165;
	t135 = 0.1e1 / t142;
	t233 = t121 * t135 * t136;
	t184 = t226 * t233;
	t195 = qJD(3) * t157 * t163;
	t241 = (-t184 + (-t153 * t166 * t190 + t158 * t195) * t136) / t131 ^ 2;
	t221 = t157 * t166;
	t138 = t158 * t198 + t218;
	t181 = t164 * t189;
	t212 = qJD(4) * t166;
	t192 = t165 * t212;
	t120 = t138 * qJD(1) - t158 * t192 + t164 * t194 - t181;
	t197 = t238 * t165;
	t141 = t158 * t219 - t197;
	t154 = 0.1e1 / t157;
	t160 = 0.1e1 / t164;
	t161 = 0.1e1 / t164 ^ 2;
	t213 = qJD(4) * t165;
	t193 = t161 * t213;
	t155 = 0.1e1 / t157 ^ 2;
	t216 = qJD(3) * t158;
	t196 = t155 * t216;
	t225 = t154 * t160;
	t240 = (t154 * t193 + t160 * t196) * t141 + t120 * t225;
	t222 = t157 * t164;
	t130 = atan2(-t138, t222);
	t125 = cos(t130);
	t124 = sin(t130);
	t232 = t124 * t138;
	t119 = t125 * t222 - t232;
	t116 = 0.1e1 / t119;
	t117 = 0.1e1 / t119 ^ 2;
	t239 = 0.2e1 * t141;
	t133 = t138 ^ 2;
	t223 = t155 * t161;
	t132 = t133 * t223 + 0.1e1;
	t128 = 0.1e1 / t132;
	t177 = t157 * t213 + t164 * t216;
	t201 = t138 * t223;
	t199 = t238 * t157;
	t182 = qJD(3) * t199;
	t183 = t165 * t190;
	t217 = qJD(1) * t166;
	t122 = t165 * t189 * t158 - t183 + (t217 * t158 - t182 - t212) * t164;
	t204 = t122 * t225;
	t108 = (t177 * t201 - t204) * t128;
	t174 = -t108 * t138 + t177;
	t104 = (-t108 * t222 - t122) * t124 + t174 * t125;
	t118 = t116 * t117;
	t237 = t104 * t118;
	t156 = t154 / t153;
	t162 = t160 * t161;
	t236 = (t122 * t201 + (-t155 * t162 * t213 - t156 * t161 * t216) * t133) / t132 ^ 2;
	t235 = t117 * t141;
	t234 = t120 * t117;
	t231 = t124 * t141;
	t230 = t124 * t157;
	t229 = t125 * t138;
	t228 = t125 * t141;
	t227 = t125 * t158;
	t224 = t155 * t158;
	t220 = t161 * t165;
	t214 = qJD(4) * t164;
	t134 = t141 ^ 2;
	t114 = t117 * t134 + 0.1e1;
	t211 = 0.2e1 * (-t134 * t237 - t141 * t234) / t114 ^ 2;
	t210 = 0.2e1 * t241;
	t209 = -0.2e1 * t236;
	t208 = t118 * t239;
	t207 = t154 * t236;
	t206 = t117 * t231;
	t202 = t138 * t225;
	t200 = t160 * t224;
	t179 = t138 * t200 + t238;
	t115 = t179 * t128;
	t191 = t238 - t115;
	t188 = t116 * t211;
	t187 = t117 * t211;
	t186 = t221 * t239;
	t185 = t160 * t207;
	t140 = t158 * t197 - t219;
	t180 = t138 * t220 - t140 * t160;
	t178 = t136 * t140 * t166 - t238 * t135;
	t126 = 0.1e1 / t131;
	t123 = t142 * qJD(1) - t158 * t181 - t165 * t182 - t192;
	t112 = 0.1e1 / t114;
	t111 = t180 * t154 * t128;
	t107 = (-t124 + (t125 * t202 + t124) * t128) * t141;
	t106 = -t115 * t229 + (t191 * t230 + t227) * t164;
	t105 = t125 * t157 * t165 - t124 * t140 + (-t124 * t222 - t229) * t111;
	t103 = t179 * t209 + (t122 * t200 + t217 + (-t193 * t224 + (-0.2e1 * t156 * t158 ^ 2 - t154) * t160 * qJD(3)) * t138) * t128;
	t101 = -0.2e1 * t180 * t207 + (-t180 * t196 + (t122 * t220 - t123 * t160 + (t140 * t220 + (-0.2e1 * t162 * t165 ^ 2 - t160) * t138) * qJD(4)) * t154) * t128;
	t1 = [t240 * t128 + t185 * t239, 0, t103, t101, 0; t138 * t188 + (-t122 * t116 + (t104 * t138 + t107 * t120) * t117) * t112 + (t107 * t187 + (0.2e1 * t107 * t237 + (t120 * t128 - t120 - (-t108 * t128 * t202 + t209) * t141) * t117 * t124 + (-(-0.2e1 * t138 * t185 - t108) * t235 + (-(t108 + t204) * t141 + t240 * t138) * t117 * t128) * t125) * t112) * t141, 0, t106 * t141 * t187 + (-(-t103 * t229 + (t108 * t232 - t122 * t125) * t115) * t235 + (t104 * t208 + t234) * t106 + (-t116 * t221 - (-t115 * t230 + t124 * t199 + t227) * t235) * t213) * t112 + (t188 * t221 + ((-t116 * t215 - (t191 * qJD(3) - t108) * t206) * t158 + (t116 * t190 + (t166 * t104 - (-t103 + t217) * t231 - (t191 * t108 - qJD(3)) * t228) * t117) * t157) * t112) * t164, (t105 * t235 - t116 * t142) * t211 + (t105 * t234 + t121 * t116 + (t105 * t208 - t117 * t142) * t104 - (t165 * t216 - t157 * t214 - t101 * t138 - t111 * t122 + (-t111 * t222 - t140) * t108) * t117 * t228 - (-t123 + (-t101 * t164 - t108 * t165) * t157 - t174 * t111) * t206) * t112, 0; t178 * t157 * t210 + (-t178 * t216 + ((qJD(1) * t135 + 0.2e1 * t140 * t233) * t166 + (-t238 * t121 - t123 * t166 + t140 * t190) * t136) * t157) * t126, 0, (t135 * t158 * t166 + t165 * t203) * t210 + (0.2e1 * t165 * t184 + t176 * t135 + ((t121 * t166 - 0.2e1 * t165 * t195) * t158 + (t163 * t214 + 0.2e1 * t166 * t183) * t153) * t136) * t126, t136 * t186 * t241 + (t186 * t233 + (t120 * t221 + (t157 * t190 - t158 * t215) * t141) * t136) * t126, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end