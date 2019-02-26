% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:27
% EndTime: 2019-02-26 20:46:28
% DurationCPUTime: 0.68s
% Computational Cost: add. (2081->92), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->90)
t139 = sin(qJ(1));
t141 = cos(qJ(1));
t135 = pkin(9) + qJ(3);
t133 = sin(t135);
t156 = qJD(1) * t133 + qJD(5);
t134 = cos(t135);
t176 = qJD(3) * t134;
t201 = t156 * t139 - t141 * t176;
t182 = t139 * t134;
t123 = atan2(-t182, t133);
t122 = cos(t123);
t121 = sin(t123);
t169 = t121 * t182;
t107 = t122 * t133 - t169;
t104 = 0.1e1 / t107;
t140 = cos(qJ(5));
t181 = t139 * t140;
t138 = sin(qJ(5));
t183 = t138 * t141;
t118 = t133 * t183 + t181;
t114 = 0.1e1 / t118;
t128 = 0.1e1 / t133;
t105 = 0.1e1 / t107 ^ 2;
t115 = 0.1e1 / t118 ^ 2;
t129 = 0.1e1 / t133 ^ 2;
t136 = t139 ^ 2;
t132 = t134 ^ 2;
t187 = t129 * t132;
t126 = t136 * t187 + 0.1e1;
t124 = 0.1e1 / t126;
t200 = t124 - 0.1e1;
t137 = t141 ^ 2;
t186 = t132 * t137;
t101 = t105 * t186 + 0.1e1;
t99 = 0.1e1 / t101;
t199 = t105 * t99;
t175 = qJD(3) * t139;
t166 = t129 * t175;
t178 = qJD(1) * t141;
t167 = t134 * t178;
t98 = ((t133 * t175 - t167) * t128 + t132 * t166) * t124;
t158 = -t98 + t175;
t159 = -t139 * t98 + qJD(3);
t190 = t122 * t134;
t93 = t159 * t190 + (t158 * t133 - t167) * t121;
t198 = t104 * t105 * t93;
t157 = qJD(5) * t133 + qJD(1);
t152 = t157 * t141;
t102 = t138 * t152 + t140 * t201;
t180 = t140 * t141;
t184 = t138 * t139;
t117 = -t133 * t180 + t184;
t113 = t117 ^ 2;
t112 = t113 * t115 + 0.1e1;
t193 = t115 * t117;
t103 = -t138 * t201 + t140 * t152;
t195 = t103 * t114 * t115;
t197 = 0.1e1 / t112 ^ 2 * (t102 * t193 - t113 * t195);
t192 = t117 * t138;
t191 = t121 * t133;
t189 = t128 * t132;
t188 = t128 * t134;
t179 = qJD(1) * t139;
t177 = qJD(3) * t133;
t154 = t132 * t139 * t178;
t174 = 0.2e1 / t101 ^ 2 * (-t186 * t198 + (-t133 * t137 * t176 - t154) * t105);
t173 = 0.2e1 * t198;
t172 = 0.2e1 * t197;
t131 = t134 * t132;
t150 = qJD(3) * (-t128 * t129 * t131 - t188);
t171 = 0.2e1 * (t129 * t154 + t136 * t150) / t126 ^ 2;
t170 = t99 * t177;
t168 = t124 * t189;
t164 = 0.1e1 + t187;
t163 = t104 * t174;
t162 = 0.2e1 * t117 * t195;
t161 = t134 * t171;
t160 = t139 * t171;
t155 = t139 * t168;
t153 = t164 * t141;
t151 = t114 * t140 + t115 * t192;
t149 = t151 * t141;
t120 = -t133 * t184 + t180;
t119 = t133 * t181 + t183;
t110 = 0.1e1 / t112;
t109 = t164 * t139 * t124;
t97 = (t200 * t134 * t121 + t122 * t155) * t141;
t95 = t139 * t191 + t190 + (-t122 * t182 - t191) * t109;
t94 = -t164 * t160 + (qJD(1) * t153 + 0.2e1 * t139 * t150) * t124;
t1 = [t128 * t141 * t161 + (qJD(3) * t153 + t179 * t188) * t124, 0, t94, 0, 0, 0; (t104 * t170 + (t163 + (qJD(1) * t97 + t93) * t199) * t134) * t139 + ((t97 * t170 + (t97 * t174 + ((t98 * t155 + t200 * t177 + t161) * t121 + (t160 * t189 + t134 * t98 + (t131 * t166 - (t98 - 0.2e1 * t175) * t134) * t124) * t122) * t99 * t141) * t134) * t105 + (t97 * t173 + (-t104 + ((t136 - t137) * t122 * t168 + t200 * t169) * t105) * qJD(1)) * t134 * t99) * t141, 0 (t104 * t99 * t179 + (t163 + (qJD(3) * t95 + t93) * t199) * t141) * t133 + (((-qJD(3) * t104 + t95 * t173) * t141 + (t95 * t179 + (-(-t109 * t178 - t139 * t94) * t122 - ((t109 * t139 - 0.1e1) * t98 + (-t109 + t139) * qJD(3)) * t121) * t134 * t141) * t105) * t99 + (t95 * t174 - ((-t94 + t178) * t121 + (t158 * t109 - t159) * t122) * t99 * t133) * t105 * t141) * t134, 0, 0, 0; (-t114 * t119 + t120 * t193) * t172 + (t120 * t162 + (-t120 * t102 - t119 * t103 + t157 * t117 * t181 - (-t134 * t175 - t156 * t141) * t192) * t115 + (t156 * t180 + (-t157 * t138 + t140 * t176) * t139) * t114) * t110, 0, t134 * t149 * t172 + (t149 * t177 + (t151 * t179 + ((qJD(5) * t114 + t162) * t138 + (-t102 * t138 + (-qJD(5) * t117 + t103) * t140) * t115) * t141) * t134) * t110, 0, -0.2e1 * t197 + 0.2e1 * (t102 * t110 * t115 + (-t110 * t195 - t115 * t197) * t117) * t117, 0;];
JaD_rot  = t1;
