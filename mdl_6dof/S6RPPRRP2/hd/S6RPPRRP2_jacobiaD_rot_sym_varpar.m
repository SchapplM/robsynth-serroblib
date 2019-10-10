% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP2
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
%   Wie in S6RPPRRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:31
	% DurationCPUTime: 1.06s
	% Computational Cost: add. (3244->94), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->93)
	t144 = qJ(1) + pkin(9);
	t140 = sin(t144);
	t134 = t140 ^ 2;
	t143 = pkin(10) + qJ(4);
	t139 = sin(t143);
	t133 = t139 ^ 2;
	t141 = cos(t143);
	t136 = 0.1e1 / t141 ^ 2;
	t193 = t133 * t136;
	t128 = t134 * t193 + 0.1e1;
	t132 = t139 * t133;
	t135 = 0.1e1 / t141;
	t190 = t135 * t139;
	t154 = qJD(4) * (t132 * t135 * t136 + t190);
	t142 = cos(t144);
	t181 = qJD(1) * t142;
	t191 = t133 * t140;
	t159 = t181 * t191;
	t197 = (t134 * t154 + t136 * t159) / t128 ^ 2;
	t207 = -0.2e1 * t197;
	t165 = 0.1e1 + t193;
	t206 = t140 * t165;
	t146 = cos(qJ(5));
	t184 = t142 * t146;
	t145 = sin(qJ(5));
	t187 = t140 * t145;
	t122 = t141 * t184 + t187;
	t162 = qJD(5) * t141 - qJD(1);
	t180 = qJD(4) * t139;
	t205 = t162 * t145 + t146 * t180;
	t188 = t140 * t139;
	t125 = atan2(-t188, -t141);
	t124 = cos(t125);
	t123 = sin(t125);
	t173 = t123 * t188;
	t109 = -t124 * t141 - t173;
	t106 = 0.1e1 / t109;
	t116 = 0.1e1 / t122;
	t107 = 0.1e1 / t109 ^ 2;
	t117 = 0.1e1 / t122 ^ 2;
	t126 = 0.1e1 / t128;
	t204 = t126 - 0.1e1;
	t172 = t126 * t133 * t135;
	t160 = t140 * t172;
	t99 = (t204 * t139 * t123 - t124 * t160) * t142;
	t203 = t107 * t99;
	t179 = qJD(4) * t140;
	t169 = t136 * t179;
	t170 = t139 * t181;
	t178 = qJD(4) * t141;
	t100 = (-(-t140 * t178 - t170) * t135 + t133 * t169) * t126;
	t194 = t124 * t139;
	t95 = (-t100 * t140 + qJD(4)) * t194 + (-t170 + (t100 - t179) * t141) * t123;
	t202 = t106 * t107 * t95;
	t155 = t141 * t187 + t184;
	t168 = t145 * t180;
	t104 = t155 * qJD(1) - t122 * qJD(5) + t142 * t168;
	t185 = t142 * t145;
	t186 = t140 * t146;
	t121 = t141 * t185 - t186;
	t115 = t121 ^ 2;
	t114 = t115 * t117 + 0.1e1;
	t196 = t117 * t121;
	t161 = -qJD(1) * t141 + qJD(5);
	t157 = t161 * t146;
	t105 = t140 * t157 - t205 * t142;
	t200 = t105 * t116 * t117;
	t201 = 0.1e1 / t114 ^ 2 * (-t104 * t196 - t115 * t200);
	t199 = t107 * t139;
	t198 = t107 * t142;
	t195 = t123 * t141;
	t138 = t142 ^ 2;
	t192 = t133 * t138;
	t189 = t139 * t142;
	t111 = t126 * t206;
	t183 = t111 - t140;
	t182 = qJD(1) * t140;
	t103 = t107 * t192 + 0.1e1;
	t177 = 0.2e1 / t103 ^ 2 * (-t192 * t202 + (t138 * t139 * t178 - t159) * t107);
	t176 = 0.2e1 * t202;
	t175 = -0.2e1 * t201;
	t174 = t121 * t200;
	t166 = t111 * t140 - 0.1e1;
	t164 = t139 * t177;
	t163 = t135 * t207;
	t158 = t165 * t142;
	t156 = -t116 * t145 + t146 * t196;
	t120 = -t141 * t186 + t185;
	t112 = 0.1e1 / t114;
	t101 = 0.1e1 / t103;
	t97 = -t140 * t195 + t194 + (-t124 * t188 + t195) * t111;
	t96 = t206 * t207 + (qJD(1) * t158 + 0.2e1 * t140 * t154) * t126;
	t1 = [t163 * t189 + (qJD(4) * t158 - t182 * t190) * t126, 0, 0, t96, 0, 0; (t106 * t164 + (-t106 * t178 + (qJD(1) * t99 + t95) * t199) * t101) * t140 + (t164 * t203 + (-t178 * t203 + (t99 * t176 + ((-t100 * t160 + 0.2e1 * t139 * t197 - t204 * t178) * t123 + (t163 * t191 + t100 * t139 + (t132 * t169 - (t100 - 0.2e1 * t179) * t139) * t126) * t124) * t198) * t139 + (-t106 + ((-t134 + t138) * t124 * t172 + t204 * t173) * t107) * t139 * qJD(1)) * t101) * t142, 0, 0, (-t106 * t141 + t97 * t199) * t142 * t177 + ((-t106 * t182 + (-qJD(4) * t97 - t95) * t198) * t141 + ((-qJD(4) * t106 + t97 * t176) * t142 + (t97 * t182 + (-(-t111 * t181 - t140 * t96) * t124 - (-t183 * qJD(4) + t166 * t100) * t123) * t189) * t107 - ((t96 - t181) * t123 + (-t166 * qJD(4) + t183 * t100) * t124) * t141 * t198) * t139) * t101, 0, 0; 0.2e1 * (t116 * t155 + t120 * t196) * t201 + (0.2e1 * t120 * t174 + (t120 * t104 + t155 * t105 + (-t205 * t140 - t142 * t157) * t121) * t117 + (t161 * t185 + (-t162 * t146 + t168) * t140) * t116) * t112, 0, 0, t156 * t175 * t189 + (t156 * t142 * t178 + (-t156 * t182 + ((-qJD(5) * t116 - 0.2e1 * t174) * t146 + (-t104 * t146 + (-qJD(5) * t121 + t105) * t145) * t117) * t142) * t139) * t112, t175 + 0.2e1 * (-t104 * t112 * t117 + (-t112 * t200 - t117 * t201) * t121) * t121, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:49:30
	% EndTime: 2019-10-09 23:49:32
	% DurationCPUTime: 1.51s
	% Computational Cost: add. (6087->126), mult. (6168->275), div. (1114->15), fcn. (7752->9), ass. (0->116)
	t168 = pkin(10) + qJ(4);
	t166 = cos(t168);
	t172 = sin(qJ(5));
	t220 = qJ(1) + pkin(9);
	t202 = sin(t220);
	t195 = t202 * t172;
	t167 = cos(t220);
	t173 = cos(qJ(5));
	t228 = t167 * t173;
	t149 = t166 * t228 + t195;
	t143 = 0.1e1 / t149 ^ 2;
	t165 = sin(t168);
	t160 = t165 ^ 2;
	t164 = t167 ^ 2;
	t235 = t160 * t164;
	t211 = t143 * t235;
	t133 = 0.1e1 + t211;
	t192 = qJD(1) * t202;
	t188 = t166 * t192;
	t191 = t202 * qJD(5);
	t223 = qJD(4) * t173;
	t205 = t165 * t223;
	t128 = (-t188 + t191) * t173 + (-t205 + (-qJD(5) * t166 + qJD(1)) * t172) * t167;
	t142 = 0.1e1 / t149;
	t242 = t128 * t142 * t143;
	t197 = t235 * t242;
	t225 = qJD(4) * t166;
	t249 = (-t197 + (-t160 * t167 * t192 + t164 * t165 * t225) * t143) / t133 ^ 2;
	t231 = t165 * t167;
	t145 = t166 * t195 + t228;
	t187 = t172 * t191;
	t221 = qJD(5) * t173;
	t203 = t167 * t221;
	t224 = qJD(4) * t167;
	t206 = t165 * t224;
	t127 = t145 * qJD(1) - t166 * t203 + t172 * t206 - t187;
	t194 = t202 * t173;
	t229 = t167 * t172;
	t148 = t166 * t229 - t194;
	t161 = 0.1e1 / t165;
	t169 = 0.1e1 / t172;
	t170 = 0.1e1 / t172 ^ 2;
	t204 = t170 * t221;
	t162 = 0.1e1 / t165 ^ 2;
	t207 = t162 * t225;
	t234 = t161 * t169;
	t248 = (t161 * t204 + t169 * t207) * t148 + t127 * t234;
	t230 = t165 * t172;
	t138 = atan2(-t145, t230);
	t135 = cos(t138);
	t134 = sin(t138);
	t241 = t134 * t145;
	t126 = t135 * t230 - t241;
	t123 = 0.1e1 / t126;
	t124 = 0.1e1 / t126 ^ 2;
	t247 = 0.2e1 * t148;
	t140 = t145 ^ 2;
	t232 = t162 * t170;
	t139 = t140 * t232 + 0.1e1;
	t136 = 0.1e1 / t139;
	t185 = t165 * t221 + t172 * t225;
	t209 = t145 * t232;
	t196 = t202 * t165;
	t189 = qJD(4) * t196;
	t190 = t173 * t192;
	t222 = qJD(5) * t172;
	t226 = qJD(1) * t167;
	t129 = -t172 * t189 - t167 * t222 - t190 + (t172 * t226 + t173 * t191) * t166;
	t212 = t129 * t234;
	t115 = (t185 * t209 - t212) * t136;
	t181 = -t115 * t145 + t185;
	t111 = (-t115 * t230 - t129) * t134 + t181 * t135;
	t125 = t123 * t124;
	t246 = t111 * t125;
	t163 = t161 / t160;
	t171 = t169 * t170;
	t245 = (t129 * t209 + (-t162 * t171 * t221 - t163 * t170 * t225) * t140) / t139 ^ 2;
	t244 = t124 * t148;
	t243 = t127 * t124;
	t240 = t134 * t148;
	t239 = t134 * t165;
	t238 = t135 * t145;
	t237 = t135 * t148;
	t236 = t135 * t166;
	t233 = t162 * t166;
	t227 = t170 * t173;
	t141 = t148 ^ 2;
	t121 = t124 * t141 + 0.1e1;
	t219 = 0.2e1 * (-t141 * t246 - t148 * t243) / t121 ^ 2;
	t218 = 0.2e1 * t249;
	t217 = -0.2e1 * t245;
	t216 = t125 * t247;
	t215 = t161 * t245;
	t214 = t124 * t240;
	t210 = t145 * t234;
	t208 = t169 * t233;
	t201 = t123 * t219;
	t200 = t124 * t219;
	t199 = t231 * t247;
	t198 = t169 * t215;
	t184 = t145 * t208 + t202;
	t122 = t184 * t136;
	t193 = t202 - t122;
	t147 = t166 * t194 - t229;
	t186 = t145 * t227 - t147 * t169;
	t183 = t143 * t147 * t167 - t202 * t142;
	t131 = 0.1e1 / t133;
	t130 = t149 * qJD(1) - t166 * t187 - t173 * t189 - t203;
	t119 = 0.1e1 / t121;
	t118 = t186 * t161 * t136;
	t114 = (-t134 + (t135 * t210 + t134) * t136) * t148;
	t113 = -t122 * t238 + (t193 * t239 + t236) * t172;
	t112 = t135 * t165 * t173 - t134 * t147 + (-t134 * t230 - t238) * t118;
	t110 = t184 * t217 + (t129 * t208 + t226 + (-t204 * t233 + (-0.2e1 * t163 * t166 ^ 2 - t161) * t169 * qJD(4)) * t145) * t136;
	t108 = -0.2e1 * t186 * t215 + (-t186 * t207 + (t129 * t227 - t130 * t169 + (t147 * t227 + (-0.2e1 * t171 * t173 ^ 2 - t169) * t145) * qJD(5)) * t161) * t136;
	t1 = [t248 * t136 + t198 * t247, 0, 0, t110, t108, 0; t145 * t201 + (-t129 * t123 + (t111 * t145 + t114 * t127) * t124) * t119 + (t114 * t200 + (0.2e1 * t114 * t246 + (t127 * t136 - t127 - (-t115 * t136 * t210 + t217) * t148) * t124 * t134 + (-(-0.2e1 * t145 * t198 - t115) * t244 + (-(t115 + t212) * t148 + t248 * t145) * t124 * t136) * t135) * t119) * t148, 0, 0, t113 * t148 * t200 + (-(-t110 * t238 + (t115 * t241 - t129 * t135) * t122) * t244 + (t111 * t216 + t243) * t113 + (-t123 * t231 - (-t122 * t239 + t134 * t196 + t236) * t244) * t221) * t119 + (t201 * t231 + ((-t123 * t224 - (t193 * qJD(4) - t115) * t214) * t166 + (t123 * t192 + (t167 * t111 - (-t110 + t226) * t240 - (t193 * t115 - qJD(4)) * t237) * t124) * t165) * t119) * t172, (t112 * t244 - t123 * t149) * t219 + (t112 * t243 + t128 * t123 + (t112 * t216 - t124 * t149) * t111 - (t166 * t223 - t165 * t222 - t108 * t145 - t118 * t129 + (-t118 * t230 - t147) * t115) * t124 * t237 - (-t130 + (-t108 * t172 - t115 * t173) * t165 - t181 * t118) * t214) * t119, 0; t183 * t165 * t218 + (-t183 * t225 + ((qJD(1) * t142 + 0.2e1 * t147 * t242) * t167 + (-t202 * t128 - t130 * t167 + t147 * t192) * t143) * t165) * t131, 0, 0, (t142 * t166 * t167 + t173 * t211) * t218 + (0.2e1 * t173 * t197 + (t188 + t206) * t142 + ((t128 * t167 - 0.2e1 * t164 * t205) * t166 + (t164 * t222 + 0.2e1 * t167 * t190) * t160) * t143) * t131, t143 * t199 * t249 + (t199 * t242 + (t127 * t231 + (t165 * t192 - t166 * t224) * t148) * t143) * t131, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end