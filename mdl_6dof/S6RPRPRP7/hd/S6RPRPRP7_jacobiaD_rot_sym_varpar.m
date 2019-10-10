% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRP7
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
%   Wie in S6RPRPRP7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 00:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:23
	% EndTime: 2019-10-10 00:39:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:24
	% EndTime: 2019-10-10 00:39:25
	% DurationCPUTime: 0.93s
	% Computational Cost: add. (2081->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->95)
	t139 = cos(qJ(1));
	t203 = 0.2e1 * t139;
	t133 = qJ(3) + pkin(9);
	t131 = sin(t133);
	t132 = cos(t133);
	t181 = t139 * t132;
	t121 = atan2(-t181, t131);
	t119 = sin(t121);
	t120 = cos(t121);
	t105 = -t119 * t181 + t120 * t131;
	t102 = 0.1e1 / t105;
	t136 = sin(qJ(5));
	t180 = t139 * t136;
	t137 = sin(qJ(1));
	t138 = cos(qJ(5));
	t182 = t137 * t138;
	t116 = t131 * t182 + t180;
	t112 = 0.1e1 / t116;
	t126 = 0.1e1 / t131;
	t103 = 0.1e1 / t105 ^ 2;
	t113 = 0.1e1 / t116 ^ 2;
	t127 = 0.1e1 / t131 ^ 2;
	t135 = t139 ^ 2;
	t130 = t132 ^ 2;
	t186 = t127 * t130;
	t124 = t135 * t186 + 0.1e1;
	t122 = 0.1e1 / t124;
	t202 = t122 - 0.1e1;
	t134 = t137 ^ 2;
	t177 = qJD(1) * t139;
	t155 = t130 * t137 * t177;
	t175 = qJD(3) * t132;
	t185 = t130 * t134;
	t174 = qJD(3) * t139;
	t165 = t127 * t174;
	t178 = qJD(1) * t137;
	t166 = t132 * t178;
	t96 = ((t131 * t174 + t166) * t126 + t130 * t165) * t122;
	t160 = -t96 + t174;
	t161 = -t139 * t96 + qJD(3);
	t189 = t120 * t132;
	t91 = t161 * t189 + (t160 * t131 + t166) * t119;
	t199 = t102 * t103 * t91;
	t99 = t103 * t185 + 0.1e1;
	t201 = (-t185 * t199 + (-t131 * t134 * t175 + t155) * t103) / t99 ^ 2;
	t97 = 0.1e1 / t99;
	t200 = t103 * t97;
	t158 = qJD(1) * t131 + qJD(5);
	t151 = t158 * t139;
	t159 = qJD(5) * t131 + qJD(1);
	t152 = t159 * t138;
	t100 = t137 * t152 + (t137 * t175 + t151) * t136;
	t179 = t139 * t138;
	t183 = t137 * t136;
	t115 = t131 * t183 - t179;
	t111 = t115 ^ 2;
	t110 = t111 * t113 + 0.1e1;
	t192 = t113 * t115;
	t153 = t159 * t136;
	t101 = t138 * t151 + (t138 * t175 - t153) * t137;
	t196 = t101 * t112 * t113;
	t198 = 0.1e1 / t110 ^ 2 * (t100 * t192 - t111 * t196);
	t129 = t132 * t130;
	t187 = t126 * t132;
	t149 = qJD(3) * (-t126 * t127 * t129 - t187);
	t194 = (-t127 * t155 + t135 * t149) / t124 ^ 2;
	t193 = t112 * t136;
	t191 = t115 * t138;
	t190 = t119 * t131;
	t188 = t126 * t130;
	t176 = qJD(3) * t131;
	t173 = -0.2e1 * t201;
	t172 = -0.2e1 * t199;
	t171 = 0.2e1 * t198;
	t170 = t102 * t201;
	t169 = t97 * t176;
	t168 = t132 * t194;
	t167 = t122 * t188;
	t164 = 0.1e1 + t186;
	t163 = t194 * t203;
	t162 = 0.2e1 * t115 * t196;
	t157 = t139 * t167;
	t156 = t202 * t132 * t119;
	t154 = t164 * t137;
	t150 = t113 * t191 - t193;
	t148 = t150 * t137;
	t147 = t132 * t174 - t158 * t137;
	t118 = t131 * t179 - t183;
	t117 = t131 * t180 + t182;
	t108 = 0.1e1 / t110;
	t107 = t164 * t139 * t122;
	t95 = (-t120 * t157 - t156) * t137;
	t93 = t139 * t190 + t189 + (-t120 * t181 - t190) * t107;
	t92 = -t164 * t163 + (-qJD(1) * t154 + t149 * t203) * t122;
	t1 = [-0.2e1 * t137 * t126 * t168 + (-qJD(3) * t154 + t177 * t187) * t122, 0, t92, 0, 0, 0; (t102 * t169 + (0.2e1 * t170 + (qJD(1) * t95 + t91) * t200) * t132) * t139 + ((-t95 * t169 + (t95 * t173 + ((t96 * t157 + t202 * t176 + 0.2e1 * t168) * t119 + (t163 * t188 + t132 * t96 + (t129 * t165 + (-t96 + 0.2e1 * t174) * t132) * t122) * t120) * t97 * t137) * t132) * t103 + (t95 * t172 + (t102 + ((t134 - t135) * t120 * t167 - t139 * t156) * t103) * qJD(1)) * t132 * t97) * t137, 0, (t102 * t97 * t177 + (-0.2e1 * t170 + (-qJD(3) * t93 - t91) * t200) * t137) * t131 + (((qJD(3) * t102 + t93 * t172) * t137 + (t93 * t177 + ((t107 * t178 - t139 * t92) * t120 + ((t107 * t139 - 0.1e1) * t96 + (-t107 + t139) * qJD(3)) * t119) * t132 * t137) * t103) * t97 + (t93 * t173 + ((-t92 - t178) * t119 + (t160 * t107 - t161) * t120) * t97 * t131) * t103 * t137) * t132, 0, 0, 0; (-t112 * t117 + t118 * t192) * t171 + (t118 * t162 + t139 * t112 * t152 + t147 * t193 + (t139 * t115 * t153 - t118 * t100 - t117 * t101 - t147 * t191) * t113) * t108, 0, t132 * t148 * t171 + (t148 * t176 + (-t150 * t177 + ((qJD(5) * t112 + t162) * t138 + (-t100 * t138 + (qJD(5) * t115 - t101) * t136) * t113) * t137) * t132) * t108, 0, -0.2e1 * t198 + 0.2e1 * (t100 * t113 * t108 + (-t108 * t196 - t113 * t198) * t115) * t115, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 00:39:24
	% EndTime: 2019-10-10 00:39:25
	% DurationCPUTime: 0.99s
	% Computational Cost: add. (2081->92), mult. (2519->205), div. (480->12), fcn. (2968->9), ass. (0->95)
	t142 = cos(qJ(1));
	t204 = 0.2e1 * t142;
	t136 = qJ(3) + pkin(9);
	t134 = sin(t136);
	t135 = cos(t136);
	t181 = t142 * t135;
	t124 = atan2(-t181, t134);
	t122 = sin(t124);
	t123 = cos(t124);
	t108 = -t122 * t181 + t123 * t134;
	t105 = 0.1e1 / t108;
	t140 = sin(qJ(1));
	t141 = cos(qJ(5));
	t182 = t140 * t141;
	t139 = sin(qJ(5));
	t184 = t139 * t142;
	t119 = t134 * t182 + t184;
	t115 = 0.1e1 / t119;
	t129 = 0.1e1 / t134;
	t106 = 0.1e1 / t108 ^ 2;
	t116 = 0.1e1 / t119 ^ 2;
	t130 = 0.1e1 / t134 ^ 2;
	t138 = t142 ^ 2;
	t133 = t135 ^ 2;
	t187 = t130 * t133;
	t127 = t138 * t187 + 0.1e1;
	t125 = 0.1e1 / t127;
	t203 = t125 - 0.1e1;
	t137 = t140 ^ 2;
	t186 = t133 * t137;
	t102 = t106 * t186 + 0.1e1;
	t178 = qJD(1) * t142;
	t158 = t133 * t140 * t178;
	t176 = qJD(3) * t135;
	t175 = qJD(3) * t142;
	t168 = t130 * t175;
	t179 = qJD(1) * t140;
	t169 = t135 * t179;
	t99 = ((t134 * t175 + t169) * t129 + t133 * t168) * t125;
	t163 = -t99 + t175;
	t164 = -t142 * t99 + qJD(3);
	t190 = t123 * t135;
	t94 = t164 * t190 + (t163 * t134 + t169) * t122;
	t200 = t105 * t106 * t94;
	t202 = 0.1e1 / t102 ^ 2 * (-t186 * t200 + (-t134 * t137 * t176 + t158) * t106);
	t159 = t203 * t135 * t122;
	t189 = t129 * t133;
	t170 = t125 * t189;
	t160 = t142 * t170;
	t98 = (-t123 * t160 - t159) * t140;
	t201 = t106 * t98;
	t161 = qJD(1) * t134 + qJD(5);
	t154 = t161 * t142;
	t162 = qJD(5) * t134 + qJD(1);
	t155 = t162 * t141;
	t103 = t140 * t155 + (t140 * t176 + t154) * t139;
	t180 = t142 * t141;
	t183 = t140 * t139;
	t118 = t134 * t183 - t180;
	t114 = t118 ^ 2;
	t113 = t114 * t116 + 0.1e1;
	t193 = t116 * t118;
	t156 = t162 * t139;
	t104 = t141 * t154 + (t141 * t176 - t156) * t140;
	t198 = t104 * t115 * t116;
	t199 = 0.1e1 / t113 ^ 2 * (t103 * t193 - t114 * t198);
	t197 = t106 * t135;
	t196 = t106 * t140;
	t132 = t135 * t133;
	t188 = t129 * t135;
	t152 = qJD(3) * (-t129 * t130 * t132 - t188);
	t195 = (-t130 * t158 + t138 * t152) / t127 ^ 2;
	t194 = t115 * t139;
	t192 = t118 * t141;
	t191 = t122 * t134;
	t177 = qJD(3) * t134;
	t174 = -0.2e1 * t200;
	t173 = 0.2e1 * t199;
	t172 = t135 * t202;
	t171 = t135 * t195;
	t167 = 0.1e1 + t187;
	t166 = 0.2e1 * t118 * t198;
	t165 = t195 * t204;
	t157 = t167 * t140;
	t153 = t116 * t192 - t194;
	t151 = t153 * t140;
	t150 = t135 * t175 - t161 * t140;
	t121 = t134 * t180 - t183;
	t120 = t134 * t184 + t182;
	t111 = 0.1e1 / t113;
	t110 = t167 * t142 * t125;
	t100 = 0.1e1 / t102;
	t96 = t142 * t191 + t190 + (-t123 * t181 - t191) * t110;
	t95 = -t167 * t165 + (-qJD(1) * t157 + t152 * t204) * t125;
	t1 = [-0.2e1 * t140 * t129 * t171 + (-qJD(3) * t157 + t178 * t188) * t125, 0, t95, 0, 0, 0; (0.2e1 * t105 * t172 + (t105 * t177 + (qJD(1) * t98 + t94) * t197) * t100) * t142 + (-0.2e1 * t172 * t201 + (-t177 * t201 + (t98 * t174 + ((t99 * t160 + t203 * t177 + 0.2e1 * t171) * t122 + (t165 * t189 + t99 * t135 + (t132 * t168 + (-t99 + 0.2e1 * t175) * t135) * t125) * t123) * t196) * t135 + (t105 + ((t137 - t138) * t123 * t170 - t142 * t159) * t106) * t135 * qJD(1)) * t100) * t140, 0, 0.2e1 * (-t105 * t134 - t96 * t197) * t140 * t202 + ((t105 * t178 + (-qJD(3) * t96 - t94) * t196) * t134 + ((qJD(3) * t105 + t96 * t174) * t140 + (t96 * t178 + ((t110 * t179 - t142 * t95) * t123 + ((t110 * t142 - 0.1e1) * t99 + (-t110 + t142) * qJD(3)) * t122) * t135 * t140) * t106 + ((-t95 - t179) * t122 + (t163 * t110 - t164) * t123) * t134 * t196) * t135) * t100, 0, 0, 0; (-t115 * t120 + t121 * t193) * t173 + (t121 * t166 + t142 * t115 * t155 + t150 * t194 + (t142 * t118 * t156 - t121 * t103 - t120 * t104 - t150 * t192) * t116) * t111, 0, t135 * t151 * t173 + (t151 * t177 + (-t153 * t178 + ((qJD(5) * t115 + t166) * t141 + (-t103 * t141 + (qJD(5) * t118 - t104) * t139) * t116) * t140) * t135) * t111, 0, -0.2e1 * t199 + 0.2e1 * (t103 * t111 * t116 + (-t111 * t198 - t116 * t199) * t118) * t118, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end