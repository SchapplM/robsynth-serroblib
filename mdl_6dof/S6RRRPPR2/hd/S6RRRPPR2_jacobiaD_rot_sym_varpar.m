% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR2
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
%   Wie in S6RRRPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:49
	% DurationCPUTime: 0.85s
	% Computational Cost: add. (4795->70), mult. (2858->159), div. (686->14), fcn. (3330->7), ass. (0->76)
	t128 = sin(qJ(1));
	t122 = t128 ^ 2;
	t120 = qJ(2) + qJ(3) + pkin(10);
	t118 = sin(t120);
	t113 = t118 ^ 2;
	t119 = cos(t120);
	t116 = 0.1e1 / t119 ^ 2;
	t166 = t113 * t116;
	t108 = t122 * t166 + 0.1e1;
	t112 = t118 * t113;
	t114 = t119 ^ 2;
	t115 = 0.1e1 / t119;
	t121 = qJD(2) + qJD(3);
	t164 = t115 * t118;
	t136 = t121 * (t112 * t115 / t114 + t164);
	t129 = cos(qJ(1));
	t156 = qJD(1) * t129;
	t148 = t128 * t156;
	t170 = 0.1e1 / t108 ^ 2 * (t122 * t136 + t148 * t166);
	t177 = -0.2e1 * t170;
	t161 = t121 * t128;
	t106 = 0.1e1 / t108;
	t149 = t118 * t156;
	t150 = t116 * t161;
	t92 = (-(-t119 * t161 - t149) * t115 + t113 * t150) * t106;
	t142 = t92 - t161;
	t160 = 0.1e1 / t128 * t129;
	t146 = 0.1e1 + t166;
	t176 = t128 * t146;
	t127 = t129 ^ 2;
	t175 = qJD(1) * (0.1e1 / t122 * t127 + 0.1e1) * t160;
	t158 = t128 * t118;
	t105 = atan2(-t158, -t119);
	t104 = cos(t105);
	t103 = sin(t105);
	t151 = t103 * t158;
	t100 = -t104 * t119 - t151;
	t97 = 0.1e1 / t100;
	t98 = 0.1e1 / t100 ^ 2;
	t174 = t106 - 0.1e1;
	t162 = t119 * t121;
	t140 = t118 * t127 * t162;
	t143 = -t128 * t92 + t121;
	t167 = t104 * t118;
	t87 = t143 * t167 + (t142 * t119 - t149) * t103;
	t172 = t87 * t97 * t98;
	t95 = t113 * t127 * t98 + 0.1e1;
	t173 = (t98 * t140 + (-t127 * t172 - t98 * t148) * t113) / t95 ^ 2;
	t93 = 0.1e1 / t95;
	t171 = t93 * t98;
	t169 = t129 * t98;
	t168 = t103 * t128;
	t165 = t113 * t128;
	t163 = t118 * t129;
	t124 = 0.1e1 / t128 ^ 2;
	t159 = t124 * t127;
	t157 = qJD(1) * t128;
	t155 = 0.2e1 * t172;
	t154 = t97 * t173;
	t111 = t114 * t159 + 0.1e1;
	t153 = 0.2e1 * (-t114 * t175 - t124 * t140) / t111 ^ 2;
	t152 = t93 * t162;
	t147 = 0.2e1 * t98 * t173;
	t145 = 0.1e1 + t159;
	t144 = t115 * t177;
	t141 = t104 * t106 * t113 * t115;
	t139 = t146 * t129;
	t138 = t145 * t118;
	t109 = 0.1e1 / t111;
	t101 = t106 * t176;
	t91 = (t174 * t118 * t103 - t128 * t141) * t129;
	t90 = t118 * t153 * t160 + (qJD(1) * t138 - t160 * t162) * t109;
	t89 = -t119 * t168 + t167 + (t103 * t119 - t104 * t158) * t101;
	t88 = t176 * t177 + (qJD(1) * t139 + 0.2e1 * t128 * t136) * t106;
	t85 = (-t97 * t93 * t157 + (-0.2e1 * t154 + (-t121 * t89 - t87) * t171) * t129) * t119 + (t89 * t129 * t147 + (-t129 * t121 * t97 - (-t104 * t128 * t88 - t142 * t103 + (-t103 * t121 - t104 * t156 + t168 * t92) * t101) * t98 * t163 + (t129 * t155 + t98 * t157) * t89 - ((t88 - t156) * t103 + (t142 * t101 + t143) * t104) * t119 * t169) * t93) * t118;
	t1 = [t144 * t163 + (t121 * t139 - t157 * t164) * t106, t88, t88, 0, 0, 0; (-t97 * t152 + (0.2e1 * t154 + (qJD(1) * t91 + t87) * t171) * t118) * t128 + (-t91 * t98 * t152 + (t91 * t147 + (t91 * t155 + ((0.2e1 * t118 * t170 + t162 + (-t115 * t92 * t165 - t162) * t106) * t103 + (t144 * t165 + t92 * t118 + (t112 * t150 - (t92 - 0.2e1 * t161) * t118) * t106) * t104) * t169) * t93) * t118 + (-t97 + (-(t122 - t127) * t141 + t174 * t151) * t98) * t118 * t93 * qJD(1)) * t129, t85, t85, 0, 0, 0; t145 * t119 * t153 + (0.2e1 * t119 * t175 + t121 * t138) * t109, t90, t90, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:18:48
	% EndTime: 2019-10-10 11:18:49
	% DurationCPUTime: 1.04s
	% Computational Cost: add. (5302->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->95)
	t180 = sin(qJ(1));
	t175 = qJ(2) + qJ(3) + pkin(10);
	t173 = sin(t175);
	t169 = 0.1e1 / t173 ^ 2;
	t174 = cos(t175);
	t172 = t174 ^ 2;
	t226 = t169 * t172;
	t201 = 0.1e1 + t226;
	t241 = t180 * t201;
	t177 = t180 ^ 2;
	t166 = t177 * t226 + 0.1e1;
	t164 = 0.1e1 / t166;
	t168 = 0.1e1 / t173;
	t182 = cos(qJ(1));
	t213 = qJD(1) * t182;
	t202 = t174 * t213;
	t176 = qJD(2) + qJD(3);
	t221 = t176 * t180;
	t204 = t169 * t221;
	t138 = ((t173 * t221 - t202) * t168 + t172 * t204) * t164;
	t240 = -t138 + t221;
	t197 = qJD(1) * t173 + qJD(6);
	t220 = t176 * t182;
	t239 = -t174 * t220 + t197 * t180;
	t219 = t180 * t174;
	t163 = atan2(-t219, t173);
	t154 = cos(t163);
	t153 = sin(t163);
	t206 = t153 * t219;
	t148 = t154 * t173 - t206;
	t145 = 0.1e1 / t148;
	t179 = sin(qJ(6));
	t216 = t182 * t179;
	t181 = cos(qJ(6));
	t217 = t180 * t181;
	t160 = t173 * t216 + t217;
	t156 = 0.1e1 / t160;
	t146 = 0.1e1 / t148 ^ 2;
	t157 = 0.1e1 / t160 ^ 2;
	t238 = t164 - 0.1e1;
	t230 = t154 * t174;
	t133 = (-t138 * t180 + t176) * t230 + (t240 * t173 - t202) * t153;
	t237 = t133 * t145 * t146;
	t198 = qJD(6) * t173 + qJD(1);
	t193 = t198 * t182;
	t143 = t179 * t193 + t239 * t181;
	t215 = t182 * t181;
	t218 = t180 * t179;
	t159 = -t173 * t215 + t218;
	t155 = t159 ^ 2;
	t152 = t155 * t157 + 0.1e1;
	t229 = t157 * t159;
	t144 = -t239 * t179 + t181 * t193;
	t234 = t144 * t156 * t157;
	t236 = (t143 * t229 - t155 * t234) / t152 ^ 2;
	t171 = t174 * t172;
	t227 = t168 * t174;
	t191 = t176 * (-t168 * t169 * t171 - t227);
	t224 = t172 * t180;
	t195 = t213 * t224;
	t235 = (t169 * t195 + t177 * t191) / t166 ^ 2;
	t233 = t146 * t174;
	t232 = t146 * t182;
	t231 = t153 * t180;
	t228 = t159 * t179;
	t178 = t182 ^ 2;
	t225 = t172 * t178;
	t223 = t173 * t176;
	t222 = t174 * t176;
	t214 = qJD(1) * t180;
	t141 = t146 * t225 + 0.1e1;
	t212 = 0.2e1 * (-t225 * t237 + (-t173 * t178 * t222 - t195) * t146) / t141 ^ 2;
	t211 = 0.2e1 * t237;
	t210 = 0.2e1 * t236;
	t209 = -0.2e1 * t235;
	t208 = t174 * t235;
	t207 = t174 * t232;
	t205 = t168 * t224;
	t200 = t174 * t212;
	t199 = 0.2e1 * t159 * t234;
	t196 = t154 * t164 * t168 * t172;
	t194 = t201 * t182;
	t192 = t156 * t181 + t157 * t228;
	t190 = t192 * t182;
	t162 = -t173 * t218 + t215;
	t161 = t173 * t217 + t216;
	t150 = 0.1e1 / t152;
	t149 = t164 * t241;
	t139 = 0.1e1 / t141;
	t137 = (t238 * t174 * t153 + t180 * t196) * t182;
	t135 = t173 * t231 + t230 + (-t153 * t173 - t154 * t219) * t149;
	t134 = t209 * t241 + (qJD(1) * t194 + 0.2e1 * t180 * t191) * t164;
	t131 = t174 * t190 * t210 + (t190 * t223 + (t192 * t214 + ((qJD(6) * t156 + t199) * t179 + (-t143 * t179 + (-qJD(6) * t159 + t144) * t181) * t157) * t182) * t174) * t150;
	t130 = (t135 * t233 + t145 * t173) * t182 * t212 + ((t145 * t214 + (t135 * t176 + t133) * t232) * t173 + (-t145 * t220 - (-t134 * t154 * t180 + t240 * t153 + (t138 * t231 - t153 * t176 - t154 * t213) * t149) * t207 + (t146 * t214 + t182 * t211) * t135 - ((-t134 + t213) * t153 + ((t149 * t180 - 0.1e1) * t176 + (-t149 + t180) * t138) * t154) * t173 * t232) * t174) * t139;
	t1 = [0.2e1 * t182 * t168 * t208 + (t176 * t194 + t214 * t227) * t164, t134, t134, 0, 0, 0; (t145 * t200 + (t145 * t223 + (qJD(1) * t137 + t133) * t233) * t139) * t180 + (t146 * t200 * t137 + (-((-0.2e1 * t208 + t223 + (-t138 * t205 - t223) * t164) * t153 + (t205 * t209 - t138 * t174 + (-t171 * t204 + (t138 - 0.2e1 * t221) * t174) * t164) * t154) * t207 + (t146 * t223 + t174 * t211) * t137 + (-t145 + ((t177 - t178) * t196 + t238 * t206) * t146) * t174 * qJD(1)) * t139) * t182, t130, t130, 0, 0, 0; (-t156 * t161 + t162 * t229) * t210 + (t162 * t199 + (-t162 * t143 - t161 * t144 + t198 * t159 * t217 - (-t176 * t219 - t197 * t182) * t228) * t157 + (t197 * t215 + (-t198 * t179 + t181 * t222) * t180) * t156) * t150, t131, t131, 0, 0, -0.2e1 * t236 + 0.2e1 * (t143 * t157 * t150 + (-t150 * t234 - t157 * t236) * t159) * t159;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end