% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR4
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
%   Wie in S6RPRRPR4_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR4_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR4_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:46
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:47
	% DurationCPUTime: 1.01s
	% Computational Cost: add. (5265->87), mult. (3414->191), div. (723->12), fcn. (4027->9), ass. (0->91)
	t152 = sin(qJ(1));
	t207 = 0.2e1 * t152;
	t146 = pkin(10) + qJ(3) + qJ(4);
	t145 = cos(t146);
	t151 = cos(pkin(11));
	t181 = t152 * t151;
	t150 = sin(pkin(11));
	t153 = cos(qJ(1));
	t185 = t150 * t153;
	t133 = t145 * t185 - t181;
	t127 = t133 ^ 2;
	t182 = t152 * t150;
	t184 = t151 * t153;
	t134 = t145 * t184 + t182;
	t129 = 0.1e1 / t134 ^ 2;
	t122 = t127 * t129 + 0.1e1;
	t131 = -t145 * t182 - t184;
	t144 = sin(t146);
	t147 = qJD(3) + qJD(4);
	t186 = t147 * t153;
	t171 = t144 * t186;
	t123 = t131 * qJD(1) - t150 * t171;
	t195 = t129 * t133;
	t132 = -t145 * t181 + t185;
	t124 = t132 * qJD(1) - t151 * t171;
	t128 = 0.1e1 / t134;
	t199 = t124 * t128 * t129;
	t206 = (t123 * t195 - t127 * t199) / t122 ^ 2;
	t148 = t152 ^ 2;
	t140 = t144 ^ 2;
	t142 = 0.1e1 / t145 ^ 2;
	t193 = t140 * t142;
	t138 = t148 * t193 + 0.1e1;
	t136 = 0.1e1 / t138;
	t141 = 0.1e1 / t145;
	t179 = qJD(1) * t153;
	t170 = t144 * t179;
	t187 = t147 * t152;
	t173 = t142 * t187;
	t110 = (-(-t145 * t187 - t170) * t141 + t140 * t173) * t136;
	t205 = t110 - t187;
	t183 = t152 * t144;
	t135 = atan2(-t183, -t145);
	t126 = cos(t135);
	t125 = sin(t135);
	t174 = t125 * t183;
	t118 = -t126 * t145 - t174;
	t115 = 0.1e1 / t118;
	t116 = 0.1e1 / t118 ^ 2;
	t204 = t136 - 0.1e1;
	t197 = t126 * t144;
	t105 = (-t110 * t152 + t147) * t197 + (t205 * t145 - t170) * t125;
	t203 = t105 * t115 * t116;
	t139 = t144 * t140;
	t190 = t141 * t144;
	t161 = t147 * (t139 * t141 * t142 + t190);
	t191 = t140 * t152;
	t164 = t179 * t191;
	t202 = (t142 * t164 + t148 * t161) / t138 ^ 2;
	t201 = t116 * t144;
	t200 = t116 * t153;
	t198 = t125 * t152;
	t196 = t128 * t150;
	t194 = t140 * t141;
	t149 = t153 ^ 2;
	t192 = t140 * t149;
	t189 = t144 * t153;
	t188 = t145 * t147;
	t180 = qJD(1) * t152;
	t113 = t116 * t192 + 0.1e1;
	t178 = 0.2e1 * (-t192 * t203 + (t144 * t149 * t188 - t164) * t116) / t113 ^ 2;
	t177 = 0.2e1 * t203;
	t176 = t116 * t189;
	t175 = t133 * t199;
	t172 = t147 * t183;
	t169 = 0.1e1 + t193;
	t168 = t144 * t178;
	t167 = -0.2e1 * t144 * t202;
	t166 = t202 * t207;
	t165 = t126 * t136 * t194;
	t163 = t169 * t153;
	t162 = t151 * t195 - t196;
	t120 = 0.1e1 / t122;
	t119 = t169 * t152 * t136;
	t111 = 0.1e1 / t113;
	t109 = (t204 * t144 * t125 - t152 * t165) * t153;
	t107 = -t145 * t198 + t197 + (t125 * t145 - t126 * t183) * t119;
	t106 = -t169 * t166 + (qJD(1) * t163 + t161 * t207) * t136;
	t103 = -0.2e1 * t162 * t189 * t206 + (t162 * t145 * t186 + (-0.2e1 * t175 * t184 + t180 * t196 + (t124 * t185 + (t123 * t153 - t133 * t180) * t151) * t129) * t144) * t120;
	t102 = (t107 * t201 - t115 * t145) * t153 * t178 + ((-t115 * t180 + (-t107 * t147 - t105) * t200) * t145 + (-t115 * t186 - (-t106 * t126 * t152 - t205 * t125 + (t110 * t198 - t125 * t147 - t126 * t179) * t119) * t176 + (t116 * t180 + t153 * t177) * t107 - ((t106 - t179) * t125 + ((-t119 * t152 + 0.1e1) * t147 + (t119 - t152) * t110) * t126) * t145 * t200) * t144) * t111;
	t1 = [t141 * t153 * t167 + (t147 * t163 - t180 * t190) * t136, 0, t106, t106, 0, 0; (t115 * t168 + (-t115 * t188 + (qJD(1) * t109 + t105) * t201) * t111) * t152 + (t116 * t168 * t109 + (-((t167 - t188 + (t110 * t141 * t191 + t188) * t136) * t125 + (t166 * t194 - t110 * t144 + (-t139 * t173 + (t110 - 0.2e1 * t187) * t144) * t136) * t126) * t176 + (-t116 * t188 + t144 * t177) * t109 + (-t115 + ((-t148 + t149) * t165 + t204 * t174) * t116) * t144 * qJD(1)) * t111) * t153, 0, t102, t102, 0, 0; 0.2e1 * (-t128 * t131 + t132 * t195) * t206 + ((-t133 * qJD(1) + t150 * t172) * t128 + 0.2e1 * t132 * t175 + (-t131 * t124 - (-t134 * qJD(1) + t151 * t172) * t133 - t132 * t123) * t129) * t120, 0, t103, t103, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:28:46
	% EndTime: 2019-10-10 01:28:47
	% DurationCPUTime: 1.12s
	% Computational Cost: add. (6073->97), mult. (3810->205), div. (753->12), fcn. (4455->9), ass. (0->96)
	t177 = sin(qJ(1));
	t238 = 0.2e1 * t177;
	t175 = t177 ^ 2;
	t172 = pkin(10) + qJ(3) + qJ(4);
	t168 = sin(t172);
	t164 = t168 ^ 2;
	t169 = cos(t172);
	t166 = 0.1e1 / t169 ^ 2;
	t223 = t164 * t166;
	t159 = t175 * t223 + 0.1e1;
	t157 = 0.1e1 / t159;
	t165 = 0.1e1 / t169;
	t178 = cos(qJ(1));
	t209 = qJD(1) * t178;
	t199 = t168 * t209;
	t174 = qJD(3) + qJD(4);
	t217 = t174 * t177;
	t202 = t166 * t217;
	t131 = (-(-t169 * t217 - t199) * t165 + t164 * t202) * t157;
	t237 = t131 - t217;
	t173 = pkin(11) + qJ(6);
	t171 = cos(t173);
	t211 = t178 * t171;
	t170 = sin(t173);
	t214 = t177 * t170;
	t153 = t169 * t211 + t214;
	t215 = t177 * t168;
	t156 = atan2(-t215, -t169);
	t155 = cos(t156);
	t154 = sin(t156);
	t203 = t154 * t215;
	t141 = -t155 * t169 - t203;
	t138 = 0.1e1 / t141;
	t147 = 0.1e1 / t153;
	t139 = 0.1e1 / t141 ^ 2;
	t148 = 0.1e1 / t153 ^ 2;
	t236 = t157 - 0.1e1;
	t225 = t155 * t168;
	t126 = (-t131 * t177 + t174) * t225 + (t237 * t169 - t199) * t154;
	t235 = t126 * t138 * t139;
	t188 = t169 * t214 + t211;
	t216 = t174 * t178;
	t200 = t168 * t216;
	t132 = t188 * qJD(1) - t153 * qJD(6) + t170 * t200;
	t212 = t178 * t170;
	t213 = t177 * t171;
	t152 = t169 * t212 - t213;
	t146 = t152 ^ 2;
	t145 = t146 * t148 + 0.1e1;
	t228 = t148 * t152;
	t193 = -qJD(1) * t169 + qJD(6);
	t194 = qJD(6) * t169 - qJD(1);
	t133 = -t194 * t212 + (t193 * t177 - t200) * t171;
	t233 = t133 * t147 * t148;
	t234 = (-t132 * t228 - t146 * t233) / t145 ^ 2;
	t163 = t168 * t164;
	t220 = t165 * t168;
	t187 = t174 * (t163 * t165 * t166 + t220);
	t221 = t164 * t177;
	t191 = t209 * t221;
	t232 = (t166 * t191 + t175 * t187) / t159 ^ 2;
	t231 = t139 * t168;
	t230 = t139 * t178;
	t229 = t147 * t170;
	t227 = t152 * t171;
	t226 = t154 * t177;
	t224 = t164 * t165;
	t176 = t178 ^ 2;
	t222 = t164 * t176;
	t219 = t168 * t178;
	t218 = t169 * t174;
	t210 = qJD(1) * t177;
	t136 = t139 * t222 + 0.1e1;
	t208 = 0.2e1 * (-t222 * t235 + (t168 * t176 * t218 - t191) * t139) / t136 ^ 2;
	t207 = 0.2e1 * t235;
	t206 = -0.2e1 * t234;
	t205 = t139 * t219;
	t204 = t152 * t233;
	t198 = 0.1e1 + t223;
	t197 = t168 * t208;
	t196 = -0.2e1 * t168 * t232;
	t195 = t232 * t238;
	t192 = t155 * t157 * t224;
	t190 = t198 * t178;
	t189 = t148 * t227 - t229;
	t186 = t174 * t215 + t193 * t178;
	t151 = -t169 * t213 + t212;
	t143 = 0.1e1 / t145;
	t142 = t198 * t177 * t157;
	t134 = 0.1e1 / t136;
	t130 = (t236 * t168 * t154 - t177 * t192) * t178;
	t129 = -t169 * t226 + t225 + (t154 * t169 - t155 * t215) * t142;
	t127 = -t198 * t195 + (qJD(1) * t190 + t187 * t238) * t157;
	t124 = t189 * t206 * t219 + (t189 * t169 * t216 + (-t189 * t210 + ((-qJD(6) * t147 - 0.2e1 * t204) * t171 + (-t132 * t171 + (-qJD(6) * t152 + t133) * t170) * t148) * t178) * t168) * t143;
	t123 = (t129 * t231 - t138 * t169) * t178 * t208 + ((-t138 * t210 + (-t129 * t174 - t126) * t230) * t169 + (-t138 * t216 - (-t127 * t155 * t177 - t237 * t154 + (t131 * t226 - t154 * t174 - t155 * t209) * t142) * t205 + (t139 * t210 + t178 * t207) * t129 - ((t127 - t209) * t154 + ((-t142 * t177 + 0.1e1) * t174 + (t142 - t177) * t131) * t155) * t169 * t230) * t168) * t134;
	t1 = [t178 * t165 * t196 + (t174 * t190 - t210 * t220) * t157, 0, t127, t127, 0, 0; (t138 * t197 + (-t138 * t218 + (qJD(1) * t130 + t126) * t231) * t134) * t177 + (t139 * t197 * t130 + (-((t196 - t218 + (t131 * t165 * t221 + t218) * t157) * t154 + (t195 * t224 - t131 * t168 + (-t163 * t202 + (t131 - 0.2e1 * t217) * t168) * t157) * t155) * t205 + (-t139 * t218 + t168 * t207) * t130 + (-t138 + ((-t175 + t176) * t192 + t236 * t203) * t139) * t168 * qJD(1)) * t134) * t178, 0, t123, t123, 0, 0; 0.2e1 * (t147 * t188 + t151 * t228) * t234 + (0.2e1 * t151 * t204 - t194 * t147 * t213 + t186 * t229 + (-t194 * t152 * t214 + t151 * t132 + t133 * t188 - t186 * t227) * t148) * t143, 0, t124, t124, 0, t206 + 0.2e1 * (-t132 * t148 * t143 + (-t143 * t233 - t148 * t234) * t152) * t152;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end