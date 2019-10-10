% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP2
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
%   Wie in S6RPRPRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:30
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:48
	% DurationCPUTime: 1.11s
	% Computational Cost: add. (3244->96), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->95)
	t146 = qJ(3) + pkin(10);
	t142 = sin(t146);
	t136 = t142 ^ 2;
	t144 = cos(t146);
	t139 = 0.1e1 / t144 ^ 2;
	t195 = t136 * t139;
	t147 = qJ(1) + pkin(9);
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
	% StartTime: 2019-10-10 00:30:47
	% EndTime: 2019-10-10 00:30:49
	% DurationCPUTime: 1.49s
	% Computational Cost: add. (6087->126), mult. (6168->275), div. (1114->15), fcn. (7752->9), ass. (0->116)
	t171 = qJ(3) + pkin(10);
	t169 = cos(t171);
	t175 = sin(qJ(5));
	t223 = qJ(1) + pkin(9);
	t205 = sin(t223);
	t198 = t205 * t175;
	t170 = cos(t223);
	t176 = cos(qJ(5));
	t231 = t170 * t176;
	t152 = t169 * t231 + t198;
	t146 = 0.1e1 / t152 ^ 2;
	t168 = sin(t171);
	t163 = t168 ^ 2;
	t167 = t170 ^ 2;
	t238 = t163 * t167;
	t214 = t146 * t238;
	t136 = 0.1e1 + t214;
	t195 = qJD(1) * t205;
	t191 = t169 * t195;
	t194 = t205 * qJD(5);
	t226 = qJD(3) * t176;
	t208 = t168 * t226;
	t131 = (-t191 + t194) * t176 + (-t208 + (-qJD(5) * t169 + qJD(1)) * t175) * t170;
	t145 = 0.1e1 / t152;
	t245 = t131 * t145 * t146;
	t200 = t238 * t245;
	t228 = qJD(3) * t169;
	t252 = (-t200 + (-t163 * t170 * t195 + t167 * t168 * t228) * t146) / t136 ^ 2;
	t234 = t168 * t170;
	t148 = t169 * t198 + t231;
	t190 = t175 * t194;
	t224 = qJD(5) * t176;
	t206 = t170 * t224;
	t227 = qJD(3) * t170;
	t209 = t168 * t227;
	t130 = t148 * qJD(1) - t169 * t206 + t175 * t209 - t190;
	t197 = t205 * t176;
	t232 = t170 * t175;
	t151 = t169 * t232 - t197;
	t164 = 0.1e1 / t168;
	t172 = 0.1e1 / t175;
	t173 = 0.1e1 / t175 ^ 2;
	t207 = t173 * t224;
	t165 = 0.1e1 / t168 ^ 2;
	t210 = t165 * t228;
	t237 = t164 * t172;
	t251 = (t164 * t207 + t172 * t210) * t151 + t130 * t237;
	t233 = t168 * t175;
	t141 = atan2(-t148, t233);
	t138 = cos(t141);
	t137 = sin(t141);
	t244 = t137 * t148;
	t129 = t138 * t233 - t244;
	t126 = 0.1e1 / t129;
	t127 = 0.1e1 / t129 ^ 2;
	t250 = 0.2e1 * t151;
	t143 = t148 ^ 2;
	t235 = t165 * t173;
	t142 = t143 * t235 + 0.1e1;
	t139 = 0.1e1 / t142;
	t188 = t168 * t224 + t175 * t228;
	t212 = t148 * t235;
	t199 = t205 * t168;
	t192 = qJD(3) * t199;
	t193 = t176 * t195;
	t225 = qJD(5) * t175;
	t229 = qJD(1) * t170;
	t132 = -t175 * t192 - t170 * t225 - t193 + (t175 * t229 + t176 * t194) * t169;
	t215 = t132 * t237;
	t118 = (t188 * t212 - t215) * t139;
	t184 = -t118 * t148 + t188;
	t114 = (-t118 * t233 - t132) * t137 + t184 * t138;
	t128 = t126 * t127;
	t249 = t114 * t128;
	t166 = t164 / t163;
	t174 = t172 * t173;
	t248 = (t132 * t212 + (-t165 * t174 * t224 - t166 * t173 * t228) * t143) / t142 ^ 2;
	t247 = t127 * t151;
	t246 = t130 * t127;
	t243 = t137 * t151;
	t242 = t137 * t168;
	t241 = t138 * t148;
	t240 = t138 * t151;
	t239 = t138 * t169;
	t236 = t165 * t169;
	t230 = t173 * t176;
	t144 = t151 ^ 2;
	t124 = t127 * t144 + 0.1e1;
	t222 = 0.2e1 * (-t144 * t249 - t151 * t246) / t124 ^ 2;
	t221 = 0.2e1 * t252;
	t220 = -0.2e1 * t248;
	t219 = t128 * t250;
	t218 = t164 * t248;
	t217 = t127 * t243;
	t213 = t148 * t237;
	t211 = t172 * t236;
	t204 = t126 * t222;
	t203 = t127 * t222;
	t202 = t234 * t250;
	t201 = t172 * t218;
	t187 = t148 * t211 + t205;
	t125 = t187 * t139;
	t196 = t205 - t125;
	t150 = t169 * t197 - t232;
	t189 = t148 * t230 - t150 * t172;
	t186 = t146 * t150 * t170 - t205 * t145;
	t134 = 0.1e1 / t136;
	t133 = t152 * qJD(1) - t169 * t190 - t176 * t192 - t206;
	t122 = 0.1e1 / t124;
	t121 = t189 * t164 * t139;
	t117 = (-t137 + (t138 * t213 + t137) * t139) * t151;
	t116 = -t125 * t241 + (t196 * t242 + t239) * t175;
	t115 = t138 * t168 * t176 - t137 * t150 + (-t137 * t233 - t241) * t121;
	t113 = t187 * t220 + (t132 * t211 + t229 + (-t207 * t236 + (-0.2e1 * t166 * t169 ^ 2 - t164) * t172 * qJD(3)) * t148) * t139;
	t111 = -0.2e1 * t189 * t218 + (-t189 * t210 + (t132 * t230 - t133 * t172 + (t150 * t230 + (-0.2e1 * t174 * t176 ^ 2 - t172) * t148) * qJD(5)) * t164) * t139;
	t1 = [t251 * t139 + t201 * t250, 0, t113, 0, t111, 0; t148 * t204 + (-t132 * t126 + (t114 * t148 + t117 * t130) * t127) * t122 + (t117 * t203 + (0.2e1 * t117 * t249 + (t130 * t139 - t130 - (-t118 * t139 * t213 + t220) * t151) * t127 * t137 + (-(-0.2e1 * t148 * t201 - t118) * t247 + (-(t118 + t215) * t151 + t251 * t148) * t127 * t139) * t138) * t122) * t151, 0, t116 * t151 * t203 + (-(-t113 * t241 + (t118 * t244 - t132 * t138) * t125) * t247 + (t114 * t219 + t246) * t116 + (-t126 * t234 - (-t125 * t242 + t137 * t199 + t239) * t247) * t224) * t122 + (t204 * t234 + ((-t126 * t227 - (t196 * qJD(3) - t118) * t217) * t169 + (t126 * t195 + (t170 * t114 - (-t113 + t229) * t243 - (t196 * t118 - qJD(3)) * t240) * t127) * t168) * t122) * t175, 0, (t115 * t247 - t126 * t152) * t222 + (t115 * t246 + t131 * t126 + (t115 * t219 - t127 * t152) * t114 - (t169 * t226 - t168 * t225 - t111 * t148 - t121 * t132 + (-t121 * t233 - t150) * t118) * t127 * t240 - (-t133 + (-t111 * t175 - t118 * t176) * t168 - t184 * t121) * t217) * t122, 0; t186 * t168 * t221 + (-t186 * t228 + ((qJD(1) * t145 + 0.2e1 * t150 * t245) * t170 + (-t205 * t131 - t133 * t170 + t150 * t195) * t146) * t168) * t134, 0, (t145 * t169 * t170 + t176 * t214) * t221 + (0.2e1 * t176 * t200 + (t191 + t209) * t145 + ((t131 * t170 - 0.2e1 * t167 * t208) * t169 + (t167 * t225 + 0.2e1 * t170 * t193) * t163) * t146) * t134, 0, t146 * t202 * t252 + (t202 * t245 + (t130 * t234 + (t168 * t195 - t169 * t227) * t151) * t146) * t134, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end