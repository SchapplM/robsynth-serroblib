% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RPPRRP7
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
%   Wie in S6RPPRRP7_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP7_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPPRRP7_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_jacobiaD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:57:59
	% EndTime: 2019-10-09 23:57:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:57:59
	% EndTime: 2019-10-09 23:57:59
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:57:59
	% EndTime: 2019-10-09 23:57:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:57:59
	% EndTime: 2019-10-09 23:57:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:57:59
	% EndTime: 2019-10-09 23:57:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:58:00
	% EndTime: 2019-10-09 23:58:01
	% DurationCPUTime: 1.01s
	% Computational Cost: add. (2081->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->95)
	t136 = cos(qJ(1));
	t200 = 0.2e1 * t136;
	t130 = pkin(9) + qJ(4);
	t128 = sin(t130);
	t129 = cos(t130);
	t177 = t136 * t129;
	t118 = atan2(-t177, t128);
	t116 = sin(t118);
	t117 = cos(t118);
	t102 = -t116 * t177 + t117 * t128;
	t99 = 0.1e1 / t102;
	t134 = sin(qJ(1));
	t135 = cos(qJ(5));
	t178 = t134 * t135;
	t133 = sin(qJ(5));
	t180 = t133 * t136;
	t113 = t128 * t178 + t180;
	t109 = 0.1e1 / t113;
	t123 = 0.1e1 / t128;
	t100 = 0.1e1 / t102 ^ 2;
	t110 = 0.1e1 / t113 ^ 2;
	t124 = 0.1e1 / t128 ^ 2;
	t132 = t136 ^ 2;
	t127 = t129 ^ 2;
	t183 = t124 * t127;
	t121 = t132 * t183 + 0.1e1;
	t119 = 0.1e1 / t121;
	t199 = t119 - 0.1e1;
	t131 = t134 ^ 2;
	t174 = qJD(1) * t136;
	t152 = t127 * t134 * t174;
	t172 = qJD(4) * t129;
	t182 = t127 * t131;
	t171 = qJD(4) * t136;
	t162 = t124 * t171;
	t175 = qJD(1) * t134;
	t163 = t129 * t175;
	t93 = ((t128 * t171 + t163) * t123 + t127 * t162) * t119;
	t157 = -t93 + t171;
	t158 = -t136 * t93 + qJD(4);
	t186 = t117 * t129;
	t88 = t158 * t186 + (t157 * t128 + t163) * t116;
	t196 = t99 * t100 * t88;
	t96 = t100 * t182 + 0.1e1;
	t198 = (-t182 * t196 + (-t128 * t131 * t172 + t152) * t100) / t96 ^ 2;
	t94 = 0.1e1 / t96;
	t197 = t100 * t94;
	t176 = t136 * t135;
	t179 = t134 * t133;
	t112 = t128 * t179 - t176;
	t108 = t112 ^ 2;
	t107 = t108 * t110 + 0.1e1;
	t189 = t110 * t112;
	t155 = qJD(1) * t128 + qJD(5);
	t148 = t155 * t136;
	t156 = qJD(5) * t128 + qJD(1);
	t150 = t156 * t133;
	t98 = t135 * t148 + (t135 * t172 - t150) * t134;
	t194 = t109 * t110 * t98;
	t149 = t156 * t135;
	t97 = t134 * t149 + (t134 * t172 + t148) * t133;
	t195 = 0.1e1 / t107 ^ 2 * (-t108 * t194 + t97 * t189);
	t126 = t129 * t127;
	t184 = t123 * t129;
	t146 = qJD(4) * (-t123 * t124 * t126 - t184);
	t191 = (-t124 * t152 + t132 * t146) / t121 ^ 2;
	t190 = t109 * t133;
	t188 = t112 * t135;
	t187 = t116 * t128;
	t185 = t123 * t127;
	t173 = qJD(4) * t128;
	t170 = -0.2e1 * t198;
	t169 = -0.2e1 * t196;
	t168 = 0.2e1 * t195;
	t167 = t99 * t198;
	t166 = t94 * t173;
	t165 = t129 * t191;
	t164 = t119 * t185;
	t161 = 0.1e1 + t183;
	t160 = 0.2e1 * t112 * t194;
	t159 = t191 * t200;
	t154 = t136 * t164;
	t153 = t199 * t129 * t116;
	t151 = t161 * t134;
	t147 = t110 * t188 - t190;
	t145 = t147 * t134;
	t144 = t129 * t171 - t155 * t134;
	t115 = t128 * t176 - t179;
	t114 = t128 * t180 + t178;
	t105 = 0.1e1 / t107;
	t104 = t161 * t136 * t119;
	t92 = (-t117 * t154 - t153) * t134;
	t90 = t136 * t187 + t186 + (-t117 * t177 - t187) * t104;
	t89 = -t161 * t159 + (-qJD(1) * t151 + t146 * t200) * t119;
	t1 = [-0.2e1 * t134 * t123 * t165 + (-qJD(4) * t151 + t174 * t184) * t119, 0, 0, t89, 0, 0; (t99 * t166 + (0.2e1 * t167 + (qJD(1) * t92 + t88) * t197) * t129) * t136 + ((-t92 * t166 + (t92 * t170 + ((t93 * t154 + t199 * t173 + 0.2e1 * t165) * t116 + (t159 * t185 + t129 * t93 + (t126 * t162 + (-t93 + 0.2e1 * t171) * t129) * t119) * t117) * t94 * t134) * t129) * t100 + (t92 * t169 + (t99 + ((t131 - t132) * t117 * t164 - t136 * t153) * t100) * qJD(1)) * t129 * t94) * t134, 0, 0, (t99 * t94 * t174 + (-0.2e1 * t167 + (-qJD(4) * t90 - t88) * t197) * t134) * t128 + (((qJD(4) * t99 + t90 * t169) * t134 + (t90 * t174 + ((t104 * t175 - t136 * t89) * t117 + ((t104 * t136 - 0.1e1) * t93 + (-t104 + t136) * qJD(4)) * t116) * t129 * t134) * t100) * t94 + (t90 * t170 + ((-t89 - t175) * t116 + (t157 * t104 - t158) * t117) * t94 * t128) * t100 * t134) * t129, 0, 0; (-t109 * t114 + t115 * t189) * t168 + (t115 * t160 + t136 * t109 * t149 + t144 * t190 + (t136 * t112 * t150 - t114 * t98 - t115 * t97 - t144 * t188) * t110) * t105, 0, 0, t129 * t145 * t168 + (t145 * t173 + (-t147 * t174 + ((qJD(5) * t109 + t160) * t135 + (-t135 * t97 + (qJD(5) * t112 - t98) * t133) * t110) * t134) * t129) * t105, -0.2e1 * t195 + 0.2e1 * (t105 * t110 * t97 + (-t105 * t194 - t110 * t195) * t112) * t112, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:58:00
	% EndTime: 2019-10-09 23:58:01
	% DurationCPUTime: 1.00s
	% Computational Cost: add. (2081->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->95)
	t139 = cos(qJ(1));
	t203 = 0.2e1 * t139;
	t133 = pkin(9) + qJ(4);
	t131 = sin(t133);
	t132 = cos(t133);
	t180 = t139 * t132;
	t121 = atan2(-t180, t131);
	t119 = sin(t121);
	t120 = cos(t121);
	t105 = -t119 * t180 + t120 * t131;
	t102 = 0.1e1 / t105;
	t137 = sin(qJ(1));
	t138 = cos(qJ(5));
	t181 = t137 * t138;
	t136 = sin(qJ(5));
	t183 = t136 * t139;
	t116 = t131 * t181 + t183;
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
	t175 = qJD(4) * t132;
	t185 = t130 * t134;
	t174 = qJD(4) * t139;
	t165 = t127 * t174;
	t178 = qJD(1) * t137;
	t166 = t132 * t178;
	t96 = ((t131 * t174 + t166) * t126 + t130 * t165) * t122;
	t160 = -t96 + t174;
	t161 = -t139 * t96 + qJD(4);
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
	t182 = t137 * t136;
	t115 = t131 * t182 - t179;
	t111 = t115 ^ 2;
	t110 = t111 * t113 + 0.1e1;
	t192 = t113 * t115;
	t153 = t159 * t136;
	t101 = t138 * t151 + (t138 * t175 - t153) * t137;
	t196 = t101 * t112 * t113;
	t198 = 0.1e1 / t110 ^ 2 * (t100 * t192 - t111 * t196);
	t129 = t132 * t130;
	t187 = t126 * t132;
	t149 = qJD(4) * (-t126 * t127 * t129 - t187);
	t194 = (-t127 * t155 + t135 * t149) / t124 ^ 2;
	t193 = t112 * t136;
	t191 = t115 * t138;
	t190 = t119 * t131;
	t188 = t126 * t130;
	t176 = qJD(4) * t131;
	t173 = -0.2e1 * t201;
	t172 = -0.2e1 * t199;
	t171 = 0.2e1 * t198;
	t170 = t102 * t201;
	t169 = t97 * t176;
	t168 = t132 * t194;
	t167 = t122 * t188;
	t164 = 0.1e1 + t186;
	t163 = 0.2e1 * t115 * t196;
	t162 = t194 * t203;
	t157 = t139 * t167;
	t156 = t202 * t132 * t119;
	t154 = t164 * t137;
	t150 = t113 * t191 - t193;
	t148 = t150 * t137;
	t147 = t132 * t174 - t158 * t137;
	t118 = t131 * t179 - t182;
	t117 = t131 * t183 + t181;
	t108 = 0.1e1 / t110;
	t107 = t164 * t139 * t122;
	t95 = (-t120 * t157 - t156) * t137;
	t93 = t139 * t190 + t189 + (-t120 * t180 - t190) * t107;
	t92 = -t164 * t162 + (-qJD(1) * t154 + t149 * t203) * t122;
	t1 = [-0.2e1 * t137 * t126 * t168 + (-qJD(4) * t154 + t177 * t187) * t122, 0, 0, t92, 0, 0; (t102 * t169 + (0.2e1 * t170 + (qJD(1) * t95 + t91) * t200) * t132) * t139 + ((-t95 * t169 + (t95 * t173 + ((t96 * t157 + t202 * t176 + 0.2e1 * t168) * t119 + (t162 * t188 + t96 * t132 + (t129 * t165 + (-t96 + 0.2e1 * t174) * t132) * t122) * t120) * t97 * t137) * t132) * t103 + (t95 * t172 + (t102 + ((t134 - t135) * t120 * t167 - t139 * t156) * t103) * qJD(1)) * t132 * t97) * t137, 0, 0, (t102 * t97 * t177 + (-0.2e1 * t170 + (-qJD(4) * t93 - t91) * t200) * t137) * t131 + (((qJD(4) * t102 + t93 * t172) * t137 + (t93 * t177 + ((t107 * t178 - t139 * t92) * t120 + ((t107 * t139 - 0.1e1) * t96 + (-t107 + t139) * qJD(4)) * t119) * t132 * t137) * t103) * t97 + (t93 * t173 + ((-t92 - t178) * t119 + (t160 * t107 - t161) * t120) * t97 * t131) * t103 * t137) * t132, 0, 0; (-t112 * t117 + t118 * t192) * t171 + (t118 * t163 + t139 * t112 * t152 + t147 * t193 + (t139 * t115 * t153 - t118 * t100 - t117 * t101 - t147 * t191) * t113) * t108, 0, 0, t132 * t148 * t171 + (t148 * t176 + (-t150 * t177 + ((qJD(5) * t112 + t163) * t138 + (-t100 * t138 + (qJD(5) * t115 - t101) * t136) * t113) * t137) * t132) * t108, -0.2e1 * t198 + 0.2e1 * (t100 * t108 * t113 + (-t108 * t196 - t113 * t198) * t115) * t115, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end