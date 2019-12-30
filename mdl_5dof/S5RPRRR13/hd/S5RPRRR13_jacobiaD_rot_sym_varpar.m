% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRRR13
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
%   Wie in S5RPRRR13_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 18:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRRR13_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRR13_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:02:48
	% EndTime: 2019-12-29 18:02:48
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:02:53
	% EndTime: 2019-12-29 18:02:53
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:02:43
	% EndTime: 2019-12-29 18:02:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:02:53
	% EndTime: 2019-12-29 18:02:53
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:02:49
	% EndTime: 2019-12-29 18:02:50
	% DurationCPUTime: 1.43s
	% Computational Cost: add. (813->90), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->96)
	t132 = cos(qJ(1));
	t195 = 0.2e1 * t132;
	t167 = qJD(3) * t132;
	t126 = t132 ^ 2;
	t128 = sin(qJ(3));
	t121 = 0.1e1 / t128 ^ 2;
	t131 = cos(qJ(3));
	t125 = t131 ^ 2;
	t179 = t121 * t125;
	t118 = t126 * t179 + 0.1e1;
	t116 = 0.1e1 / t118;
	t120 = 0.1e1 / t128;
	t158 = t121 * t167;
	t129 = sin(qJ(1));
	t171 = qJD(1) * t131;
	t160 = t129 * t171;
	t90 = ((t128 * t167 + t160) * t120 + t125 * t158) * t116;
	t152 = -t90 + t167;
	t173 = t132 * t131;
	t115 = atan2(-t173, t128);
	t113 = sin(t115);
	t114 = cos(t115);
	t99 = -t113 * t173 + t114 * t128;
	t96 = 0.1e1 / t99;
	t127 = sin(qJ(4));
	t175 = t132 * t127;
	t130 = cos(qJ(4));
	t177 = t129 * t130;
	t110 = t128 * t177 + t175;
	t106 = 0.1e1 / t110;
	t107 = 0.1e1 / t110 ^ 2;
	t97 = 0.1e1 / t99 ^ 2;
	t194 = t116 - 0.1e1;
	t123 = t129 ^ 2;
	t170 = qJD(1) * t132;
	t147 = t125 * t129 * t170;
	t169 = qJD(3) * t128;
	t163 = t97 * t169;
	t153 = -t132 * t90 + qJD(3);
	t181 = t114 * t131;
	t85 = t153 * t181 + (t152 * t128 + t160) * t113;
	t192 = t85 * t96 * t97;
	t95 = t123 * t125 * t97 + 0.1e1;
	t193 = (t97 * t147 + (-t125 * t192 - t131 * t163) * t123) / t95 ^ 2;
	t93 = 0.1e1 / t95;
	t191 = t93 * t97;
	t190 = t96 * t93;
	t174 = t132 * t130;
	t178 = t129 * t127;
	t109 = t128 * t178 - t174;
	t105 = t109 ^ 2;
	t104 = t105 * t107 + 0.1e1;
	t184 = t107 * t109;
	t150 = qJD(1) * t128 + qJD(4);
	t143 = t150 * t132;
	t151 = qJD(4) * t128 + qJD(1);
	t145 = t151 * t127;
	t168 = qJD(3) * t131;
	t92 = t130 * t143 + (t130 * t168 - t145) * t129;
	t188 = t106 * t107 * t92;
	t144 = t151 * t130;
	t91 = t129 * t144 + (t129 * t168 + t143) * t127;
	t189 = 0.1e1 / t104 ^ 2 * (-t105 * t188 + t91 * t184);
	t187 = t129 * t97;
	t124 = t131 * t125;
	t141 = qJD(3) * (-t121 * t124 - t131) * t120;
	t186 = (-t121 * t147 + t126 * t141) / t118 ^ 2;
	t185 = t106 * t127;
	t183 = t109 * t130;
	t182 = t113 * t132;
	t180 = t120 * t125;
	t176 = t129 * t131;
	t172 = qJD(1) * t129;
	t166 = -0.2e1 * t192;
	t165 = 0.2e1 * t189;
	t164 = t96 * t193;
	t162 = t131 * t186;
	t161 = t116 * t180;
	t159 = t131 * t170;
	t157 = -0.2e1 * t97 * t193;
	t156 = 0.1e1 + t179;
	t155 = 0.2e1 * t109 * t188;
	t154 = t186 * t195;
	t149 = t132 * t161;
	t148 = t194 * t131 * t113;
	t146 = t156 * t129;
	t142 = t107 * t183 - t185;
	t140 = -t150 * t129 + t131 * t167;
	t112 = t128 * t174 - t178;
	t111 = t128 * t175 + t177;
	t103 = t156 * t132 * t116;
	t101 = 0.1e1 / t104;
	t89 = (-t114 * t149 - t148) * t129;
	t88 = t128 * t182 + t181 + (-t113 * t128 - t114 * t173) * t103;
	t86 = -t156 * t154 + (-qJD(1) * t146 + t141 * t195) * t116;
	t1 = [-0.2e1 * t129 * t120 * t162 + (-qJD(3) * t146 + t120 * t159) * t116, 0, t86, 0, 0; (t169 * t190 + (0.2e1 * t164 + (qJD(1) * t89 + t85) * t191) * t131) * t132 + (t89 * t157 * t131 + (-t89 * t163 + (t89 * t166 + ((t90 * t149 + t194 * t169 + 0.2e1 * t162) * t113 + (t154 * t180 + t131 * t90 + (t124 * t158 + (-t90 + 0.2e1 * t167) * t131) * t116) * t114) * t187) * t131 + (t96 + ((t123 - t126) * t114 * t161 - t132 * t148) * t97) * t171) * t93) * t129, 0, (t170 * t190 + (-0.2e1 * t164 + (-qJD(3) * t88 - t85) * t191) * t129) * t128 + (t88 * t129 * t157 + (t129 * qJD(3) * t96 + (-t114 * t132 * t86 + t152 * t113 + (-qJD(3) * t113 + t114 * t172 + t182 * t90) * t103) * t97 * t176 + (t129 * t166 + t97 * t170) * t88 + ((-t86 - t172) * t113 + (t152 * t103 - t153) * t114) * t128 * t187) * t93) * t131, 0, 0; (-t106 * t111 + t112 * t184) * t165 + (t112 * t155 + t132 * t106 * t144 + t140 * t185 + (t132 * t109 * t145 - t111 * t92 - t112 * t91 - t140 * t183) * t107) * t101, 0, t142 * t165 * t176 + (-t142 * t159 + (t142 * t169 + ((qJD(4) * t106 + t155) * t130 + (-t130 * t91 + (qJD(4) * t109 - t92) * t127) * t107) * t131) * t129) * t101, -0.2e1 * t189 + 0.2e1 * (t91 * t107 * t101 + (-t101 * t188 - t107 * t189) * t109) * t109, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 18:02:48
	% EndTime: 2019-12-29 18:02:50
	% DurationCPUTime: 1.59s
	% Computational Cost: add. (1381->94), mult. (2734->205), div. (498->12), fcn. (3199->9), ass. (0->97)
	t150 = sin(qJ(3));
	t143 = 0.1e1 / t150 ^ 2;
	t152 = cos(qJ(3));
	t147 = t152 ^ 2;
	t197 = t143 * t147;
	t153 = cos(qJ(1));
	t215 = 0.2e1 * t153;
	t214 = t152 * t197;
	t190 = t153 * t152;
	t134 = atan2(-t190, t150);
	t132 = sin(t134);
	t133 = cos(t134);
	t121 = -t132 * t190 + t133 * t150;
	t118 = 0.1e1 / t121;
	t149 = qJ(4) + qJ(5);
	t139 = sin(t149);
	t140 = cos(t149);
	t151 = sin(qJ(1));
	t193 = t151 * t140;
	t129 = t139 * t153 + t150 * t193;
	t125 = 0.1e1 / t129;
	t142 = 0.1e1 / t150;
	t119 = 0.1e1 / t121 ^ 2;
	t126 = 0.1e1 / t129 ^ 2;
	t148 = t153 ^ 2;
	t137 = t148 * t197 + 0.1e1;
	t135 = 0.1e1 / t137;
	t213 = t135 - 0.1e1;
	t145 = t151 ^ 2;
	t196 = t145 * t147;
	t114 = t119 * t196 + 0.1e1;
	t187 = qJD(1) * t153;
	t168 = t147 * t151 * t187;
	t185 = qJD(3) * t152;
	t188 = qJD(1) * t152;
	t177 = t151 * t188;
	t184 = qJD(3) * t153;
	t111 = ((t150 * t184 + t177) * t142 + t184 * t197) * t135;
	t199 = t133 * t152;
	t105 = (-t111 * t153 + qJD(3)) * t199 + (t177 + (-t111 + t184) * t150) * t132;
	t210 = t105 * t118 * t119;
	t212 = (-t196 * t210 + (-t145 * t150 * t185 + t168) * t119) / t114 ^ 2;
	t141 = qJD(4) + qJD(5);
	t171 = qJD(1) * t150 + t141;
	t161 = t151 * t185 + t171 * t153;
	t172 = t141 * t150 + qJD(1);
	t165 = t140 * t172;
	t109 = t161 * t139 + t151 * t165;
	t191 = t153 * t140;
	t194 = t151 * t139;
	t128 = t150 * t194 - t191;
	t124 = t128 ^ 2;
	t117 = t124 * t126 + 0.1e1;
	t201 = t126 * t128;
	t166 = t139 * t172;
	t110 = t161 * t140 - t151 * t166;
	t209 = t110 * t125 * t126;
	t211 = (t109 * t201 - t124 * t209) / t117 ^ 2;
	t208 = t111 * t132;
	t207 = t111 * t152;
	t206 = t119 * t151;
	t205 = t119 * t152;
	t163 = qJD(3) * (-t152 - t214) * t142;
	t204 = (-t143 * t168 + t148 * t163) / t137 ^ 2;
	t175 = 0.1e1 + t197;
	t123 = t175 * t153 * t135;
	t203 = t123 * t153;
	t202 = t125 * t139;
	t200 = t128 * t140;
	t198 = t142 * t147;
	t195 = t150 * t153;
	t192 = t151 * t152;
	t189 = qJD(1) * t151;
	t186 = qJD(3) * t150;
	t183 = 0.2e1 * t211;
	t182 = -0.2e1 * t210;
	t181 = t152 * t212;
	t180 = t119 * t192;
	t179 = t152 * t204;
	t178 = t135 * t198;
	t176 = t152 * t187;
	t174 = 0.2e1 * t128 * t209;
	t173 = t204 * t215;
	t170 = t153 * t178;
	t169 = t213 * t152 * t132;
	t167 = t175 * t151;
	t164 = t126 * t200 - t202;
	t162 = -t171 * t151 + t152 * t184;
	t131 = t150 * t191 - t194;
	t130 = t139 * t195 + t193;
	t115 = 0.1e1 / t117;
	t112 = 0.1e1 / t114;
	t108 = (-t133 * t170 - t169) * t151;
	t107 = t132 * t195 + t199 + (-t132 * t150 - t133 * t190) * t123;
	t106 = -t175 * t173 + (-qJD(1) * t167 + t163 * t215) * t135;
	t102 = -0.2e1 * t211 + 0.2e1 * (t109 * t115 * t126 + (-t115 * t209 - t126 * t211) * t128) * t128;
	t1 = [-0.2e1 * t151 * t142 * t179 + (-qJD(3) * t167 + t142 * t176) * t135, 0, t106, 0, 0; (0.2e1 * t118 * t181 + (t118 * t186 + (qJD(1) * t108 + t105) * t205) * t112) * t153 + (-0.2e1 * t119 * t181 * t108 + (((t111 * t170 + t213 * t186 + 0.2e1 * t179) * t132 + (t173 * t198 + t207 + (-t207 + (0.2e1 * t152 + t214) * t184) * t135) * t133) * t180 + (-t119 * t186 + t152 * t182) * t108 + (t118 + ((t145 - t148) * t133 * t178 - t153 * t169) * t119) * t188) * t112) * t151, 0, 0.2e1 * (-t107 * t205 - t118 * t150) * t151 * t212 + ((t118 * t187 + (-qJD(3) * t107 - t105) * t206) * t150 + (t151 * qJD(3) * t118 + (-t106 * t133 * t153 + t132 * t184 + t203 * t208 - t208 + (-qJD(3) * t132 + t133 * t189) * t123) * t180 + (t119 * t187 + t151 * t182) * t107 + ((-t106 - t189) * t132 + ((-0.1e1 + t203) * qJD(3) + (-t123 + t153) * t111) * t133) * t150 * t206) * t152) * t112, 0, 0; (-t125 * t130 + t131 * t201) * t183 + (t131 * t174 + t153 * t125 * t165 + t162 * t202 + (t153 * t128 * t166 - t131 * t109 - t130 * t110 - t162 * t200) * t126) * t115, 0, t164 * t183 * t192 + (-t164 * t176 + (t164 * t186 + ((t125 * t141 + t174) * t140 + (-t109 * t140 + (t128 * t141 - t110) * t139) * t126) * t152) * t151) * t115, t102, t102;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end