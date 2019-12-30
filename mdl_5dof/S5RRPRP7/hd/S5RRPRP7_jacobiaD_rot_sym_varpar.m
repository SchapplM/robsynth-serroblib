% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP7
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
%   Wie in S5RRPRP7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPRP7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:50:30
	% EndTime: 2019-12-29 18:50:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:50:30
	% EndTime: 2019-12-29 18:50:30
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:50:30
	% EndTime: 2019-12-29 18:50:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:50:35
	% EndTime: 2019-12-29 18:50:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:50:30
	% EndTime: 2019-12-29 18:50:32
	% DurationCPUTime: 1.53s
	% Computational Cost: add. (2270->94), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->93)
	t139 = sin(qJ(1));
	t136 = t139 ^ 2;
	t135 = qJ(2) + pkin(8);
	t133 = sin(t135);
	t129 = t133 ^ 2;
	t134 = cos(t135);
	t131 = 0.1e1 / t134 ^ 2;
	t188 = t129 * t131;
	t124 = t136 * t188 + 0.1e1;
	t128 = t133 * t129;
	t130 = 0.1e1 / t134;
	t185 = t130 * t133;
	t149 = qJD(2) * (t128 * t130 * t131 + t185);
	t141 = cos(qJ(1));
	t177 = qJD(1) * t141;
	t186 = t129 * t139;
	t154 = t177 * t186;
	t194 = (t131 * t154 + t136 * t149) / t124 ^ 2;
	t204 = -0.2e1 * t194;
	t161 = 0.1e1 + t188;
	t203 = t139 * t161;
	t140 = cos(qJ(4));
	t179 = t141 * t140;
	t138 = sin(qJ(4));
	t182 = t139 * t138;
	t120 = t134 * t179 + t182;
	t183 = t139 * t133;
	t121 = atan2(-t183, -t134);
	t116 = cos(t121);
	t115 = sin(t121);
	t168 = t115 * t183;
	t105 = -t116 * t134 - t168;
	t102 = 0.1e1 / t105;
	t112 = 0.1e1 / t120;
	t103 = 0.1e1 / t105 ^ 2;
	t113 = 0.1e1 / t120 ^ 2;
	t122 = 0.1e1 / t124;
	t202 = t122 - 0.1e1;
	t137 = t141 ^ 2;
	t176 = qJD(2) * t134;
	t187 = t129 * t137;
	t175 = qJD(2) * t139;
	t163 = t131 * t175;
	t164 = t133 * t177;
	t96 = (-(-t134 * t175 - t164) * t130 + t129 * t163) * t122;
	t158 = t96 - t175;
	t159 = -t139 * t96 + qJD(2);
	t190 = t116 * t133;
	t91 = t159 * t190 + (t158 * t134 - t164) * t115;
	t199 = t102 * t103 * t91;
	t99 = t103 * t187 + 0.1e1;
	t201 = (-t187 * t199 + (t133 * t137 * t176 - t154) * t103) / t99 ^ 2;
	t97 = 0.1e1 / t99;
	t200 = t103 * t97;
	t150 = t134 * t182 + t179;
	t174 = qJD(2) * t141;
	t162 = t133 * t174;
	t100 = t150 * qJD(1) - t120 * qJD(4) + t138 * t162;
	t180 = t141 * t138;
	t181 = t139 * t140;
	t119 = t134 * t180 - t181;
	t111 = t119 ^ 2;
	t110 = t111 * t113 + 0.1e1;
	t192 = t113 * t119;
	t156 = -qJD(1) * t134 + qJD(4);
	t157 = qJD(4) * t134 - qJD(1);
	t101 = -t157 * t180 + (t156 * t139 - t162) * t140;
	t196 = t101 * t112 * t113;
	t198 = 0.1e1 / t110 ^ 2 * (-t100 * t192 - t111 * t196);
	t193 = t112 * t138;
	t191 = t115 * t134;
	t189 = t119 * t140;
	t184 = t133 * t141;
	t178 = qJD(1) * t139;
	t173 = 0.2e1 * t201;
	t172 = 0.2e1 * t199;
	t171 = -0.2e1 * t198;
	t170 = t102 * t201;
	t169 = t97 * t176;
	t167 = t119 * t196;
	t166 = t122 * t129 * t130;
	t160 = t130 * t204;
	t155 = t139 * t166;
	t153 = t161 * t141;
	t152 = t156 * t141;
	t151 = t113 * t189 - t193;
	t118 = -t134 * t181 + t180;
	t108 = 0.1e1 / t110;
	t107 = t122 * t203;
	t95 = (t202 * t133 * t115 - t116 * t155) * t141;
	t93 = -t139 * t191 + t190 + (-t116 * t183 + t191) * t107;
	t92 = t203 * t204 + (qJD(1) * t153 + 0.2e1 * t139 * t149) * t122;
	t1 = [t160 * t184 + (qJD(2) * t153 - t178 * t185) * t122, t92, 0, 0, 0; (-t102 * t169 + (0.2e1 * t170 + (qJD(1) * t95 + t91) * t200) * t133) * t139 + ((-t95 * t169 + (t95 * t173 + ((0.2e1 * t133 * t194 - t96 * t155 - t202 * t176) * t115 + (t160 * t186 + t133 * t96 + (t128 * t163 - (t96 - 0.2e1 * t175) * t133) * t122) * t116) * t97 * t141) * t133) * t103 + (t95 * t172 + (-t102 + ((-t136 + t137) * t116 * t166 + t202 * t168) * t103) * qJD(1)) * t133 * t97) * t141, (-t102 * t97 * t178 + (-0.2e1 * t170 + (-qJD(2) * t93 - t91) * t200) * t141) * t134 + (((-qJD(2) * t102 + t93 * t172) * t141 + (t93 * t178 + (-(-t107 * t177 - t139 * t92) * t116 - ((t107 * t139 - 0.1e1) * t96 + (-t107 + t139) * qJD(2)) * t115) * t184) * t103) * t97 + (t93 * t173 - ((t92 - t177) * t115 + (t158 * t107 + t159) * t116) * t97 * t134) * t103 * t141) * t133, 0, 0, 0; 0.2e1 * (t112 * t150 + t118 * t192) * t198 + (0.2e1 * t118 * t167 - t157 * t112 * t181 + (t133 * t175 + t152) * t193 + (t118 * t100 + t150 * t101 - t152 * t189 - (qJD(2) * t133 * t140 + t157 * t138) * t119 * t139) * t113) * t108, t151 * t171 * t184 + (t151 * t134 * t174 + (-t151 * t178 + ((-qJD(4) * t112 - 0.2e1 * t167) * t140 + (-t100 * t140 + (-qJD(4) * t119 + t101) * t138) * t113) * t141) * t133) * t108, 0, t171 + 0.2e1 * (-t100 * t113 * t108 + (-t108 * t196 - t113 * t198) * t119) * t119, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:50:30
	% EndTime: 2019-12-29 18:50:33
	% DurationCPUTime: 2.22s
	% Computational Cost: add. (3877->124), mult. (6168->272), div. (1114->15), fcn. (7752->9), ass. (0->115)
	t162 = qJ(2) + pkin(8);
	t161 = cos(t162);
	t167 = sin(qJ(4));
	t241 = sin(qJ(1));
	t201 = t241 * t167;
	t168 = cos(qJ(4));
	t169 = cos(qJ(1));
	t221 = t169 * t168;
	t145 = t161 * t221 + t201;
	t139 = 0.1e1 / t145 ^ 2;
	t160 = sin(t162);
	t156 = t160 ^ 2;
	t166 = t169 ^ 2;
	t229 = t156 * t166;
	t206 = t139 * t229;
	t134 = 0.1e1 + t206;
	t193 = qJD(1) * t241;
	t218 = qJD(2) * t169;
	t197 = t160 * t218;
	t179 = t161 * t193 + t197;
	t192 = t241 * qJD(4);
	t222 = t169 * t167;
	t124 = (-qJD(4) * t161 + qJD(1)) * t222 + (t192 - t179) * t168;
	t138 = 0.1e1 / t145;
	t236 = t124 * t138 * t139;
	t187 = t229 * t236;
	t198 = qJD(2) * t160 * t166;
	t244 = (-t187 + (-t156 * t169 * t193 + t161 * t198) * t139) / t134 ^ 2;
	t224 = t160 * t169;
	t141 = t161 * t201 + t221;
	t184 = t167 * t192;
	t215 = qJD(4) * t169;
	t195 = t168 * t215;
	t123 = t141 * qJD(1) - t161 * t195 + t167 * t197 - t184;
	t200 = t241 * t168;
	t144 = t161 * t222 - t200;
	t157 = 0.1e1 / t160;
	t163 = 0.1e1 / t167;
	t164 = 0.1e1 / t167 ^ 2;
	t216 = qJD(4) * t168;
	t196 = t164 * t216;
	t158 = 0.1e1 / t160 ^ 2;
	t219 = qJD(2) * t161;
	t199 = t158 * t219;
	t228 = t157 * t163;
	t243 = (t157 * t196 + t163 * t199) * t144 + t123 * t228;
	t225 = t160 * t167;
	t133 = atan2(-t141, t225);
	t128 = cos(t133);
	t127 = sin(t133);
	t235 = t127 * t141;
	t122 = t128 * t225 - t235;
	t119 = 0.1e1 / t122;
	t120 = 0.1e1 / t122 ^ 2;
	t242 = 0.2e1 * t144;
	t136 = t141 ^ 2;
	t226 = t158 * t164;
	t135 = t136 * t226 + 0.1e1;
	t131 = 0.1e1 / t135;
	t180 = t160 * t216 + t167 * t219;
	t204 = t141 * t226;
	t202 = t241 * t160;
	t185 = qJD(2) * t202;
	t186 = t168 * t193;
	t220 = qJD(1) * t169;
	t125 = t168 * t192 * t161 - t186 + (t220 * t161 - t185 - t215) * t167;
	t207 = t125 * t228;
	t111 = (t180 * t204 - t207) * t131;
	t177 = -t111 * t141 + t180;
	t107 = (-t111 * t225 - t125) * t127 + t177 * t128;
	t121 = t119 * t120;
	t240 = t107 * t121;
	t159 = t157 / t156;
	t165 = t163 * t164;
	t239 = (t125 * t204 + (-t158 * t165 * t216 - t159 * t164 * t219) * t136) / t135 ^ 2;
	t238 = t120 * t144;
	t237 = t123 * t120;
	t234 = t127 * t144;
	t233 = t127 * t160;
	t232 = t128 * t141;
	t231 = t128 * t144;
	t230 = t128 * t161;
	t227 = t158 * t161;
	t223 = t164 * t168;
	t217 = qJD(4) * t167;
	t137 = t144 ^ 2;
	t117 = t120 * t137 + 0.1e1;
	t214 = 0.2e1 * (-t137 * t240 - t144 * t237) / t117 ^ 2;
	t213 = 0.2e1 * t244;
	t212 = -0.2e1 * t239;
	t211 = t121 * t242;
	t210 = t157 * t239;
	t209 = t120 * t234;
	t205 = t141 * t228;
	t203 = t163 * t227;
	t182 = t141 * t203 + t241;
	t118 = t182 * t131;
	t194 = t241 - t118;
	t191 = t119 * t214;
	t190 = t120 * t214;
	t189 = t224 * t242;
	t188 = t163 * t210;
	t143 = t161 * t200 - t222;
	t183 = t141 * t223 - t143 * t163;
	t181 = t139 * t143 * t169 - t241 * t138;
	t129 = 0.1e1 / t134;
	t126 = t145 * qJD(1) - t161 * t184 - t168 * t185 - t195;
	t115 = 0.1e1 / t117;
	t114 = t183 * t157 * t131;
	t110 = (-t127 + (t128 * t205 + t127) * t131) * t144;
	t109 = -t118 * t232 + (t194 * t233 + t230) * t167;
	t108 = t128 * t160 * t168 - t127 * t143 + (-t127 * t225 - t232) * t114;
	t106 = t182 * t212 + (t125 * t203 + t220 + (-t196 * t227 + (-0.2e1 * t159 * t161 ^ 2 - t157) * t163 * qJD(2)) * t141) * t131;
	t104 = -0.2e1 * t183 * t210 + (-t183 * t199 + (t125 * t223 - t126 * t163 + (t143 * t223 + (-0.2e1 * t165 * t168 ^ 2 - t163) * t141) * qJD(4)) * t157) * t131;
	t1 = [t243 * t131 + t188 * t242, t106, 0, t104, 0; t141 * t191 + (-t125 * t119 + (t107 * t141 + t110 * t123) * t120) * t115 + (t110 * t190 + (0.2e1 * t110 * t240 + (t123 * t131 - t123 - (-t111 * t131 * t205 + t212) * t144) * t120 * t127 + (-(-0.2e1 * t141 * t188 - t111) * t238 + (-(t111 + t207) * t144 + t243 * t141) * t120 * t131) * t128) * t115) * t144, t109 * t144 * t190 + (-(-t106 * t232 + (t111 * t235 - t125 * t128) * t118) * t238 + (t107 * t211 + t237) * t109 + (-t119 * t224 - (-t118 * t233 + t127 * t202 + t230) * t238) * t216) * t115 + (t191 * t224 + ((-t119 * t218 - (t194 * qJD(2) - t111) * t209) * t161 + (t119 * t193 + (t169 * t107 - (-t106 + t220) * t234 - (t111 * t194 - qJD(2)) * t231) * t120) * t160) * t115) * t167, 0, (t108 * t238 - t119 * t145) * t214 + (t108 * t237 + t124 * t119 + (t108 * t211 - t120 * t145) * t107 - (t168 * t219 - t160 * t217 - t104 * t141 - t114 * t125 + (-t114 * t225 - t143) * t111) * t120 * t231 - (-t126 + (-t104 * t167 - t111 * t168) * t160 - t177 * t114) * t209) * t115, 0; t181 * t160 * t213 + (-t181 * t219 + ((qJD(1) * t138 + 0.2e1 * t143 * t236) * t169 + (-t241 * t124 - t126 * t169 + t143 * t193) * t139) * t160) * t129, (t138 * t161 * t169 + t168 * t206) * t213 + (0.2e1 * t168 * t187 + t179 * t138 + ((t124 * t169 - 0.2e1 * t168 * t198) * t161 + (t166 * t217 + 0.2e1 * t169 * t186) * t156) * t139) * t129, 0, t139 * t189 * t244 + (t189 * t236 + (t123 * t224 + (t160 * t193 - t161 * t218) * t144) * t139) * t129, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end