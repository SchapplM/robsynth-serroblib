% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRR3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:10
% EndTime: 2019-02-26 21:55:11
% DurationCPUTime: 0.71s
% Computational Cost: add. (2270->94), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->93)
t139 = sin(qJ(1));
t136 = t139 ^ 2;
t135 = qJ(2) + pkin(11);
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
t91 = t159 * t190 + (t134 * t158 - t164) * t115;
t199 = t102 * t103 * t91;
t99 = t103 * t187 + 0.1e1;
t201 = (-t187 * t199 + (t133 * t137 * t176 - t154) * t103) / t99 ^ 2;
t97 = 0.1e1 / t99;
t200 = t103 * t97;
t150 = t134 * t182 + t179;
t174 = qJD(2) * t141;
t162 = t133 * t174;
t100 = qJD(1) * t150 - t120 * qJD(4) + t138 * t162;
t180 = t141 * t138;
t181 = t139 * t140;
t119 = t134 * t180 - t181;
t111 = t119 ^ 2;
t110 = t111 * t113 + 0.1e1;
t192 = t113 * t119;
t156 = -qJD(1) * t134 + qJD(4);
t157 = qJD(4) * t134 - qJD(1);
t101 = -t157 * t180 + (t139 * t156 - t162) * t140;
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
t151 = t189 * t113 - t193;
t118 = -t134 * t181 + t180;
t108 = 0.1e1 / t110;
t107 = t122 * t203;
t95 = (t115 * t133 * t202 - t116 * t155) * t141;
t93 = -t139 * t191 + t190 + (-t116 * t183 + t191) * t107;
t92 = t203 * t204 + (qJD(1) * t153 + 0.2e1 * t139 * t149) * t122;
t1 = [t160 * t184 + (qJD(2) * t153 - t178 * t185) * t122, t92, 0, 0, 0, 0; (-t102 * t169 + (0.2e1 * t170 + (qJD(1) * t95 + t91) * t200) * t133) * t139 + ((-t95 * t169 + (t95 * t173 + ((0.2e1 * t133 * t194 - t155 * t96 - t176 * t202) * t115 + (t160 * t186 + t133 * t96 + (t128 * t163 - (t96 - 0.2e1 * t175) * t133) * t122) * t116) * t97 * t141) * t133) * t103 + (t95 * t172 + (-t102 + ((-t136 + t137) * t116 * t166 + t202 * t168) * t103) * qJD(1)) * t133 * t97) * t141 (-t102 * t97 * t178 + (-0.2e1 * t170 + (-qJD(2) * t93 - t91) * t200) * t141) * t134 + (((-qJD(2) * t102 + t172 * t93) * t141 + (t93 * t178 + (-(-t107 * t177 - t139 * t92) * t116 - ((t107 * t139 - 0.1e1) * t96 + (-t107 + t139) * qJD(2)) * t115) * t184) * t103) * t97 + (t93 * t173 - ((t92 - t177) * t115 + (t107 * t158 + t159) * t116) * t97 * t134) * t103 * t141) * t133, 0, 0, 0, 0; 0.2e1 * (t112 * t150 + t118 * t192) * t198 + (0.2e1 * t118 * t167 - t157 * t112 * t181 + (t133 * t175 + t152) * t193 + (t118 * t100 + t150 * t101 - t152 * t189 - (qJD(2) * t133 * t140 + t138 * t157) * t119 * t139) * t113) * t108, t151 * t171 * t184 + (t151 * t134 * t174 + (-t151 * t178 + ((-qJD(4) * t112 - 0.2e1 * t167) * t140 + (-t100 * t140 + (-qJD(4) * t119 + t101) * t138) * t113) * t141) * t133) * t108, 0, t171 + 0.2e1 * (-t100 * t113 * t108 + (-t108 * t196 - t113 * t198) * t119) * t119, 0, 0;];
JaD_rot  = t1;
