% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP9
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
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP9_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:01
% EndTime: 2019-02-26 20:48:02
% DurationCPUTime: 0.66s
% Computational Cost: add. (1161->91), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->95)
t135 = cos(qJ(1));
t198 = 0.2e1 * t135;
t124 = pkin(9) + qJ(5);
t122 = sin(t124);
t132 = sin(qJ(3));
t123 = cos(t124);
t133 = sin(qJ(1));
t179 = t133 * t123;
t112 = t122 * t135 + t132 * t179;
t109 = 0.1e1 / t112 ^ 2;
t177 = t135 * t123;
t180 = t133 * t122;
t111 = t132 * t180 - t177;
t186 = t111 * t123;
t108 = 0.1e1 / t112;
t188 = t108 * t122;
t146 = t109 * t186 - t188;
t107 = t111 ^ 2;
t100 = t107 * t109 + 0.1e1;
t97 = 0.1e1 / t100;
t197 = t146 * t97;
t134 = cos(qJ(3));
t176 = t135 * t134;
t117 = atan2(-t176, t132);
t115 = sin(t117);
t116 = cos(t117);
t104 = -t115 * t176 + t116 * t132;
t101 = 0.1e1 / t104;
t125 = 0.1e1 / t132;
t102 = 0.1e1 / t104 ^ 2;
t126 = 0.1e1 / t132 ^ 2;
t131 = t135 ^ 2;
t130 = t134 ^ 2;
t183 = t126 * t130;
t120 = t131 * t183 + 0.1e1;
t118 = 0.1e1 / t120;
t196 = t118 - 0.1e1;
t128 = t133 ^ 2;
t173 = qJD(1) * t135;
t150 = t130 * t133 * t173;
t171 = qJD(3) * t134;
t182 = t128 * t130;
t170 = qJD(3) * t135;
t160 = t126 * t170;
t174 = qJD(1) * t134;
t162 = t133 * t174;
t94 = ((t132 * t170 + t162) * t125 + t130 * t160) * t118;
t155 = -t94 + t170;
t156 = -t135 * t94 + qJD(3);
t185 = t116 * t134;
t88 = t156 * t185 + (t155 * t132 + t162) * t115;
t192 = t101 * t102 * t88;
t99 = t102 * t182 + 0.1e1;
t195 = (-t182 * t192 + (-t128 * t132 * t171 + t150) * t102) / t99 ^ 2;
t187 = t109 * t111;
t153 = qJD(1) * t132 + qJD(5);
t144 = t133 * t171 + t153 * t135;
t154 = qJD(5) * t132 + qJD(1);
t148 = t122 * t154;
t93 = t144 * t123 - t133 * t148;
t191 = t108 * t109 * t93;
t147 = t123 * t154;
t92 = t144 * t122 + t133 * t147;
t194 = (-t107 * t191 + t92 * t187) / t100 ^ 2;
t95 = 0.1e1 / t99;
t193 = t102 * t95;
t129 = t134 * t130;
t145 = qJD(3) * (-t126 * t129 - t134) * t125;
t189 = (-t126 * t150 + t131 * t145) / t120 ^ 2;
t184 = t125 * t130;
t181 = t132 * t135;
t175 = qJD(1) * t133;
t172 = qJD(3) * t132;
t169 = -0.2e1 * t195;
t168 = 0.2e1 * t194;
t167 = -0.2e1 * t192;
t166 = t101 * t195;
t165 = t95 * t172;
t164 = t134 * t189;
t163 = t118 * t184;
t161 = t134 * t173;
t159 = 0.1e1 + t183;
t158 = 0.2e1 * t111 * t191;
t157 = t189 * t198;
t152 = t135 * t163;
t151 = t196 * t134 * t115;
t149 = t159 * t133;
t143 = -t153 * t133 + t134 * t170;
t114 = t132 * t177 - t180;
t113 = t122 * t181 + t179;
t106 = t159 * t135 * t118;
t91 = (-t116 * t152 - t151) * t133;
t90 = t115 * t181 + t185 + (-t115 * t132 - t116 * t176) * t106;
t89 = -t159 * t157 + (-qJD(1) * t149 + t145 * t198) * t118;
t1 = [-0.2e1 * t133 * t125 * t164 + (-qJD(3) * t149 + t125 * t161) * t118, 0, t89, 0, 0, 0; (t101 * t165 + (0.2e1 * t166 + (qJD(1) * t91 + t88) * t193) * t134) * t135 + (t91 * t134 * t95 * t167 + (-t91 * t165 + (t91 * t169 + ((t94 * t152 + t196 * t172 + 0.2e1 * t164) * t115 + (t157 * t184 + t94 * t134 + (t129 * t160 + (-t94 + 0.2e1 * t170) * t134) * t118) * t116) * t95 * t133) * t134) * t102 + (t101 + ((t128 - t131) * t116 * t163 - t135 * t151) * t102) * t95 * t174) * t133, 0 (t101 * t95 * t173 + (-0.2e1 * t166 + (-qJD(3) * t90 - t88) * t193) * t133) * t132 + (((qJD(3) * t101 + t90 * t167) * t133 + (t90 * t173 + ((t106 * t175 - t135 * t89) * t116 + ((t106 * t135 - 0.1e1) * t94 + (-t106 + t135) * qJD(3)) * t115) * t133 * t134) * t102) * t95 + (t90 * t169 + ((-t89 - t175) * t115 + (t155 * t106 - t156) * t116) * t95 * t132) * t102 * t133) * t134, 0, 0, 0; (-t108 * t113 + t114 * t187) * t168 + (t114 * t158 + t135 * t108 * t147 + t143 * t188 + (t135 * t111 * t148 - t113 * t93 - t114 * t92 - t143 * t186) * t109) * t97, 0, -t161 * t197 + (t172 * t197 + (t146 * t168 + ((qJD(5) * t108 + t158) * t123 + (-t123 * t92 + (qJD(5) * t111 - t93) * t122) * t109) * t97) * t134) * t133, 0, -0.2e1 * t194 + 0.2e1 * (t109 * t92 * t97 + (-t109 * t194 - t97 * t191) * t111) * t111, 0;];
JaD_rot  = t1;
