% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRR7
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
%   Wie in S5RPRRR7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRR7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:04:51
	% EndTime: 2019-12-31 19:04:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:04:51
	% EndTime: 2019-12-31 19:04:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:04:51
	% EndTime: 2019-12-31 19:04:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:04:51
	% EndTime: 2019-12-31 19:04:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:04:52
	% EndTime: 2019-12-31 19:04:52
	% DurationCPUTime: 0.69s
	% Computational Cost: add. (1976->91), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->94)
	t123 = qJ(1) + pkin(9);
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
	t1 = [t151 * t176 + (qJD(3) * t144 - t126 * t158) * t113, 0, t83, 0, 0; (-t167 * t189 + (0.2e1 * t164 + (qJD(1) * t86 + t82) * t190) * t130) * t121 + (t86 * t153 * t130 + (-t86 * t162 + (t86 * t166 + ((0.2e1 * t130 * t187 - t87 * t146 - t194 * t167) * t110 + (t151 * t180 + t130 * t87 + (t124 * t156 - (t87 - 0.2e1 * t169) * t130) * t113) * t111) * t185) * t130 + (-t93 + (-(t119 - t120) * t111 * t160 + t194 * t161) * t94) * t170) * t90) * t122, 0, (-t172 * t189 + (-0.2e1 * t164 + (-qJD(3) * t85 - t82) * t190) * t122) * t132 + (t85 * t122 * t153 + (-t122 * qJD(3) * t93 - ((-t100 * t171 - t121 * t83) * t111 + (-t150 * t100 - t149) * t110) * t94 * t176 + (t122 * t166 + t94 * t172) * t85 - ((t83 - t171) * t110 + (t149 * t100 + t150) * t111) * t132 * t185) * t90) * t130, 0, 0; 0.2e1 * (t103 * t141 + t107 * t183) * t191 + (0.2e1 * t107 * t163 + (t141 * t89 + t107 * t88 + (-t195 * t121 - t122 * t143) * t108) * t104 + (t147 * t177 + (-t148 * t131 + t155) * t121) * t103) * t98, 0, -t158 * t197 + (t167 * t197 + (t142 * t165 + ((-qJD(4) * t103 - 0.2e1 * t163) * t131 + (-t131 * t88 + (-qJD(4) * t108 + t89) * t129) * t104) * t98) * t130) * t122, t165 + 0.2e1 * (-t104 * t88 * t98 + (-t104 * t191 - t98 * t188) * t108) * t108, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:04:52
	% EndTime: 2019-12-31 19:04:52
	% DurationCPUTime: 0.74s
	% Computational Cost: add. (2635->97), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
	t158 = sin(qJ(3));
	t153 = t158 ^ 2;
	t159 = cos(qJ(3));
	t155 = 0.1e1 / t159 ^ 2;
	t198 = t153 * t155;
	t151 = qJ(1) + pkin(9);
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
	t1 = [t147 * t154 * t177 + (-t146 * t154 * t195 + qJD(3) * t172) * t140, 0, t109, 0, 0; (t121 * t179 + (-t121 * t192 + (qJD(1) * t111 + t108) * t211) * t115) * t146 + (t122 * t179 * t111 + (-((t114 * t174 + t218 * t192 + t177) * t136 + (t178 * t199 - t213 + (t213 + (t219 - t221) * t194) * t140) * t137) * t187 + (-t122 * t192 + t158 * t189) * t111 + (-t121 + ((-t144 + t145) * t137 * t185 + t218 * t186) * t122) * t195) * t115) * t147, 0, (t110 * t211 - t121 * t159) * t147 * t191 + ((-t121 * t197 + (-qJD(3) * t110 - t108) * t212) * t159 + (-t147 * qJD(3) * t121 - (-t109 * t137 * t146 + t136 * t194 + t209 * t214 - t214 + (-qJD(3) * t136 - t137 * t196) * t126) * t187 + (t122 * t197 + t147 * t189) * t110 - ((t109 - t196) * t136 + ((0.1e1 - t209) * qJD(3) + (t126 - t146) * t114) * t137) * t159 * t212) * t158) * t115, 0, 0; (t128 * t169 + t132 * t208) * t190 + (0.2e1 * t132 * t188 + (t132 * t112 + t169 * t113 + (-t146 * t220 - t147 * t171) * t133) * t129 + (t175 * t202 + (-t176 * t149 + t182) * t146) * t128) * t118, 0, -t147 * t167 * t190 + (-t167 * t197 + (t170 * t192 + ((-t128 * t150 - 0.2e1 * t188) * t149 + (-t112 * t149 + (-t133 * t150 + t113) * t148) * t129) * t158) * t147) * t118, t105, t105;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end