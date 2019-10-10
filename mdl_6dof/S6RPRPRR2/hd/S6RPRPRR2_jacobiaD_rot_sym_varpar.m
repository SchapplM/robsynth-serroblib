% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR2
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
%   Wie in S6RPRPRR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_jacobiaD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:59
	% DurationCPUTime: 1.10s
	% Computational Cost: add. (3244->96), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->95)
	t146 = qJ(3) + pkin(11);
	t142 = sin(t146);
	t136 = t142 ^ 2;
	t144 = cos(t146);
	t139 = 0.1e1 / t144 ^ 2;
	t195 = t136 * t139;
	t147 = qJ(1) + pkin(10);
	t143 = sin(t147);
	t213 = 0.2e1 * t143;
	t212 = t142 * t195;
	t145 = cos(t147);
	t149 = cos(qJ(5));
	t187 = t145 * t149;
	t148 = sin(qJ(5));
	t190 = t143 * t148;
	t125 = t144 * t187 + t190;
	t165 = qJD(5) * t144 - qJD(1);
	t184 = qJD(3) * t142;
	t211 = t165 * t148 + t149 * t184;
	t191 = t143 * t142;
	t128 = atan2(-t191, -t144);
	t127 = cos(t128);
	t126 = sin(t128);
	t175 = t126 * t191;
	t112 = -t127 * t144 - t175;
	t109 = 0.1e1 / t112;
	t119 = 0.1e1 / t125;
	t138 = 0.1e1 / t144;
	t110 = 0.1e1 / t112 ^ 2;
	t120 = 0.1e1 / t125 ^ 2;
	t210 = -0.2e1 * t142;
	t137 = t143 ^ 2;
	t131 = t137 * t195 + 0.1e1;
	t129 = 0.1e1 / t131;
	t209 = t129 - 0.1e1;
	t185 = qJD(1) * t145;
	t172 = t142 * t185;
	t182 = qJD(3) * t144;
	t183 = qJD(3) * t143;
	t103 = (-(-t143 * t182 - t172) * t138 + t183 * t195) * t129;
	t197 = t127 * t142;
	t98 = (-t103 * t143 + qJD(3)) * t197 + (-t172 + (t103 - t183) * t144) * t126;
	t208 = t109 * t110 * t98;
	t158 = t144 * t190 + t187;
	t171 = t148 * t184;
	t107 = t158 * qJD(1) - qJD(5) * t125 + t145 * t171;
	t188 = t145 * t148;
	t189 = t143 * t149;
	t124 = t144 * t188 - t189;
	t118 = t124 ^ 2;
	t117 = t118 * t120 + 0.1e1;
	t199 = t120 * t124;
	t164 = -qJD(1) * t144 + qJD(5);
	t160 = t164 * t149;
	t108 = t143 * t160 - t145 * t211;
	t204 = t108 * t119 * t120;
	t207 = (-t107 * t199 - t118 * t204) / t117 ^ 2;
	t206 = t103 * t126;
	t205 = t103 * t142;
	t203 = t110 * t142;
	t202 = t110 * t145;
	t193 = t138 * t142;
	t157 = qJD(3) * (t138 * t212 + t193);
	t162 = t136 * t143 * t185;
	t201 = (t137 * t157 + t139 * t162) / t131 ^ 2;
	t169 = 0.1e1 + t195;
	t114 = t169 * t143 * t129;
	t200 = t114 * t143;
	t198 = t126 * t144;
	t196 = t136 * t138;
	t141 = t145 ^ 2;
	t194 = t136 * t141;
	t192 = t142 * t145;
	t186 = qJD(1) * t143;
	t181 = qJD(3) * t145;
	t106 = t110 * t194 + 0.1e1;
	t180 = 0.2e1 / t106 ^ 2 * (-t194 * t208 + (t141 * t142 * t182 - t162) * t110);
	t179 = 0.2e1 * t208;
	t178 = -0.2e1 * t207;
	t177 = t124 * t204;
	t176 = t110 * t192;
	t174 = t129 * t196;
	t168 = t142 * t180;
	t167 = t201 * t210;
	t166 = t201 * t213;
	t163 = t143 * t174;
	t161 = t169 * t145;
	t159 = -t119 * t148 + t149 * t199;
	t123 = -t144 * t189 + t188;
	t115 = 0.1e1 / t117;
	t104 = 0.1e1 / t106;
	t102 = (t209 * t142 * t126 - t127 * t163) * t145;
	t100 = -t143 * t198 + t197 + (-t127 * t191 + t198) * t114;
	t99 = -t169 * t166 + (qJD(1) * t161 + t157 * t213) * t129;
	t1 = [t138 * t145 * t167 + (qJD(3) * t161 - t186 * t193) * t129, 0, t99, 0, 0, 0; (t109 * t168 + (-t109 * t182 + (qJD(1) * t102 + t98) * t203) * t104) * t143 + (t110 * t168 * t102 + (-((t103 * t163 + t209 * t182 + t167) * t126 + (t166 * t196 - t205 + (t205 + (t210 - t212) * t183) * t129) * t127) * t176 + (-t110 * t182 + t142 * t179) * t102 + (-t109 + ((-t137 + t141) * t127 * t174 + t209 * t175) * t110) * t142 * qJD(1)) * t104) * t145, 0, (t100 * t203 - t109 * t144) * t145 * t180 + ((-t109 * t186 + (-qJD(3) * t100 - t98) * t202) * t144 + (-t109 * t181 - (-t127 * t143 * t99 + t126 * t183 + t200 * t206 - t206 + (-qJD(3) * t126 - t127 * t185) * t114) * t176 + (t110 * t186 + t145 * t179) * t100 - ((t99 - t185) * t126 + ((0.1e1 - t200) * qJD(3) + (t114 - t143) * t103) * t127) * t144 * t202) * t142) * t104, 0, 0, 0; 0.2e1 * (t119 * t158 + t123 * t199) * t207 + (0.2e1 * t123 * t177 + (t123 * t107 + t158 * t108 + (-t143 * t211 - t145 * t160) * t124) * t120 + (t164 * t188 + (-t165 * t149 + t171) * t143) * t119) * t115, 0, t159 * t178 * t192 + (t159 * t144 * t181 + (-t159 * t186 + ((-qJD(5) * t119 - 0.2e1 * t177) * t149 + (-t107 * t149 + (-qJD(5) * t124 + t108) * t148) * t120) * t145) * t142) * t115, 0, t178 + 0.2e1 * (-t107 * t115 * t120 + (-t115 * t204 - t120 * t207) * t124) * t124, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:47:58
	% EndTime: 2019-10-10 00:47:59
	% DurationCPUTime: 1.12s
	% Computational Cost: add. (3952->98), mult. (2734->203), div. (498->12), fcn. (3199->9), ass. (0->98)
	t176 = qJ(3) + pkin(11);
	t169 = sin(t176);
	t163 = t169 ^ 2;
	t171 = cos(t176);
	t166 = 0.1e1 / t171 ^ 2;
	t224 = t163 * t166;
	t177 = qJ(1) + pkin(10);
	t170 = sin(t177);
	t242 = 0.2e1 * t170;
	t241 = t169 * t224;
	t172 = cos(t177);
	t178 = qJ(5) + qJ(6);
	t174 = cos(t178);
	t216 = t172 * t174;
	t173 = sin(t178);
	t219 = t170 * t173;
	t152 = t171 * t216 + t219;
	t175 = qJD(5) + qJD(6);
	t194 = t171 * t175 - qJD(1);
	t213 = qJD(3) * t169;
	t240 = t194 * t173 + t174 * t213;
	t220 = t170 * t169;
	t155 = atan2(-t220, -t171);
	t154 = cos(t155);
	t153 = sin(t155);
	t204 = t153 * t220;
	t139 = -t154 * t171 - t204;
	t136 = 0.1e1 / t139;
	t146 = 0.1e1 / t152;
	t165 = 0.1e1 / t171;
	t137 = 0.1e1 / t139 ^ 2;
	t147 = 0.1e1 / t152 ^ 2;
	t239 = -0.2e1 * t169;
	t164 = t170 ^ 2;
	t159 = t164 * t224 + 0.1e1;
	t157 = 0.1e1 / t159;
	t238 = t157 - 0.1e1;
	t214 = qJD(1) * t172;
	t201 = t169 * t214;
	t211 = qJD(3) * t171;
	t212 = qJD(3) * t170;
	t130 = (-(-t170 * t211 - t201) * t165 + t212 * t224) * t157;
	t226 = t154 * t169;
	t125 = (-t130 * t170 + qJD(3)) * t226 + (-t201 + (t130 - t212) * t171) * t153;
	t237 = t125 * t136 * t137;
	t187 = t171 * t219 + t216;
	t200 = t173 * t213;
	t131 = t187 * qJD(1) - t152 * t175 + t172 * t200;
	t217 = t172 * t173;
	t218 = t170 * t174;
	t151 = t171 * t217 - t218;
	t145 = t151 ^ 2;
	t144 = t145 * t147 + 0.1e1;
	t228 = t147 * t151;
	t193 = -qJD(1) * t171 + t175;
	t189 = t193 * t174;
	t132 = t170 * t189 - t240 * t172;
	t233 = t132 * t146 * t147;
	t236 = (-t131 * t228 - t145 * t233) / t144 ^ 2;
	t235 = t130 * t153;
	t234 = t130 * t169;
	t232 = t137 * t169;
	t231 = t137 * t172;
	t222 = t165 * t169;
	t186 = qJD(3) * (t165 * t241 + t222);
	t191 = t163 * t170 * t214;
	t230 = (t164 * t186 + t166 * t191) / t159 ^ 2;
	t198 = 0.1e1 + t224;
	t141 = t198 * t170 * t157;
	t229 = t141 * t170;
	t227 = t153 * t171;
	t225 = t163 * t165;
	t168 = t172 ^ 2;
	t223 = t163 * t168;
	t221 = t169 * t172;
	t215 = qJD(1) * t170;
	t210 = qJD(3) * t172;
	t135 = t137 * t223 + 0.1e1;
	t209 = 0.2e1 * (-t223 * t237 + (t168 * t169 * t211 - t191) * t137) / t135 ^ 2;
	t208 = 0.2e1 * t237;
	t207 = -0.2e1 * t236;
	t206 = t151 * t233;
	t205 = t137 * t221;
	t203 = t157 * t225;
	t197 = t169 * t209;
	t196 = t230 * t239;
	t195 = t230 * t242;
	t192 = t170 * t203;
	t190 = t198 * t172;
	t188 = -t146 * t173 + t174 * t228;
	t150 = -t171 * t218 + t217;
	t142 = 0.1e1 / t144;
	t133 = 0.1e1 / t135;
	t129 = (t238 * t169 * t153 - t154 * t192) * t172;
	t128 = -t170 * t227 + t226 + (-t154 * t220 + t227) * t141;
	t126 = -t198 * t195 + (qJD(1) * t190 + t186 * t242) * t157;
	t123 = t207 + 0.2e1 * (-t131 * t142 * t147 + (-t142 * t233 - t147 * t236) * t151) * t151;
	t1 = [t165 * t172 * t196 + (qJD(3) * t190 - t215 * t222) * t157, 0, t126, 0, 0, 0; (t136 * t197 + (-t136 * t211 + (qJD(1) * t129 + t125) * t232) * t133) * t170 + (t137 * t197 * t129 + (-((t130 * t192 + t238 * t211 + t196) * t153 + (t195 * t225 - t234 + (t234 + (t239 - t241) * t212) * t157) * t154) * t205 + (-t137 * t211 + t169 * t208) * t129 + (-t136 + ((-t164 + t168) * t154 * t203 + t238 * t204) * t137) * t169 * qJD(1)) * t133) * t172, 0, (t128 * t232 - t136 * t171) * t172 * t209 + ((-t136 * t215 + (-qJD(3) * t128 - t125) * t231) * t171 + (-t136 * t210 - (-t126 * t154 * t170 + t153 * t212 + t229 * t235 - t235 + (-qJD(3) * t153 - t154 * t214) * t141) * t205 + (t137 * t215 + t172 * t208) * t128 - ((t126 - t214) * t153 + ((0.1e1 - t229) * qJD(3) + (t141 - t170) * t130) * t154) * t171 * t231) * t169) * t133, 0, 0, 0; 0.2e1 * (t146 * t187 + t150 * t228) * t236 + (0.2e1 * t150 * t206 + (t150 * t131 + t187 * t132 + (-t170 * t240 - t172 * t189) * t151) * t147 + (t193 * t217 + (-t194 * t174 + t200) * t170) * t146) * t142, 0, t188 * t207 * t221 + (t188 * t171 * t210 + (-t188 * t215 + ((-t146 * t175 - 0.2e1 * t206) * t174 + (-t131 * t174 + (-t151 * t175 + t132) * t173) * t147) * t172) * t169) * t142, 0, t123, t123;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end