% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:57
% EndTime: 2019-02-26 20:46:57
% DurationCPUTime: 0.68s
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
JaD_rot  = t1;
