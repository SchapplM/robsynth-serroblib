% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:54
% EndTime: 2019-02-26 20:56:55
% DurationCPUTime: 0.67s
% Computational Cost: add. (1618->95), mult. (2545->213), div. (484->12), fcn. (2996->9), ass. (0->99)
t127 = qJ(1) + pkin(9);
t125 = sin(t127);
t123 = t125 ^ 2;
t134 = sin(qJ(3));
t129 = t134 ^ 2;
t136 = cos(qJ(3));
t131 = 0.1e1 / t136 ^ 2;
t180 = t129 * t131;
t120 = t123 * t180 + 0.1e1;
t128 = t134 * t129;
t130 = 0.1e1 / t136;
t144 = qJD(3) * (t128 * t131 + t134) * t130;
t126 = cos(t127);
t176 = qJD(1) * t126;
t185 = t125 * t129;
t148 = t176 * t185;
t193 = 0.1e1 / t120 ^ 2 * (t123 * t144 + t131 * t148);
t202 = -0.2e1 * t193;
t133 = sin(qJ(4));
t179 = t133 * t136;
t135 = cos(qJ(4));
t183 = t125 * t135;
t113 = -t126 * t179 + t183;
t107 = t113 ^ 2;
t178 = t135 * t136;
t114 = t125 * t133 + t126 * t178;
t109 = 0.1e1 / t114 ^ 2;
t191 = t107 * t109;
t102 = 0.1e1 + t191;
t189 = t109 * t113;
t150 = qJD(1) * t136 - qJD(4);
t146 = t133 * t150;
t171 = qJD(4) * t136;
t151 = -qJD(1) + t171;
t173 = qJD(3) * t134;
t198 = -t133 * t173 + t151 * t135;
t93 = t125 * t146 - t198 * t126;
t166 = t93 * t189;
t108 = 0.1e1 / t114;
t110 = t108 * t109;
t190 = t107 * t110;
t156 = t135 * t173;
t161 = t125 * t178;
t94 = qJD(1) * t161 - qJD(4) * t183 - t133 * t176 + (t133 * t171 + t156) * t126;
t201 = 0.1e1 / t102 ^ 2 * (t94 * t190 + t166);
t145 = t108 * t133 + t135 * t189;
t99 = 0.1e1 / t102;
t200 = t145 * t99;
t155 = 0.1e1 + t180;
t199 = t125 * t155;
t184 = t125 * t134;
t119 = atan2(t184, t136);
t116 = cos(t119);
t115 = sin(t119);
t163 = t115 * t184;
t106 = t116 * t136 + t163;
t103 = 0.1e1 / t106;
t104 = 0.1e1 / t106 ^ 2;
t117 = 0.1e1 / t120;
t197 = t117 - 0.1e1;
t124 = t126 ^ 2;
t172 = qJD(3) * t136;
t186 = t124 * t129;
t174 = qJD(3) * t125;
t158 = t131 * t174;
t175 = qJD(1) * t134;
t159 = t126 * t175;
t92 = ((t125 * t172 + t159) * t130 + t129 * t158) * t117;
t152 = -t92 + t174;
t153 = t125 * t92 - qJD(3);
t187 = t116 * t134;
t87 = t153 * t187 + (t152 * t136 + t159) * t115;
t194 = t103 * t104 * t87;
t97 = t104 * t186 + 0.1e1;
t196 = (-t186 * t194 + (t124 * t134 * t172 - t148) * t104) / t97 ^ 2;
t95 = 0.1e1 / t97;
t195 = t104 * t95;
t188 = t115 * t136;
t182 = t126 * t134;
t181 = t126 * t135;
t177 = qJD(1) * t125;
t170 = -0.2e1 * t196;
t169 = 0.2e1 * t201;
t168 = -0.2e1 * t194;
t167 = t103 * t196;
t165 = t113 * t110 * t94;
t164 = t95 * t172;
t162 = t117 * t129 * t130;
t160 = t125 * t175;
t154 = t130 * t202;
t149 = t125 * t162;
t147 = t155 * t126;
t112 = t126 * t133 - t161;
t111 = t125 * t179 + t181;
t101 = t117 * t199;
t91 = (-t197 * t134 * t115 + t116 * t149) * t126;
t90 = t125 * t188 - t187 + (t116 * t184 - t188) * t101;
t88 = t199 * t202 + (qJD(1) * t147 + 0.2e1 * t125 * t144) * t117;
t1 = [t154 * t182 + (qJD(3) * t147 - t130 * t160) * t117, 0, t88, 0, 0, 0; (t103 * t164 + (-0.2e1 * t167 + (-qJD(1) * t91 - t87) * t195) * t134) * t125 + (t91 * t134 * t95 * t168 + (t91 * t164 + (t91 * t170 + ((0.2e1 * t134 * t193 - t92 * t149 - t197 * t172) * t115 + (t154 * t185 + t134 * t92 + (t128 * t158 + (-t92 + 0.2e1 * t174) * t134) * t117) * t116) * t95 * t126) * t134) * t104 + (t103 + ((-t123 + t124) * t116 * t162 + t197 * t163) * t104) * t95 * t175) * t126, 0 (t103 * t95 * t177 + (0.2e1 * t167 + (qJD(3) * t90 + t87) * t195) * t126) * t136 + (((qJD(3) * t103 + t90 * t168) * t126 + (-t90 * t177 + ((t101 * t176 + t125 * t88) * t116 + ((-t101 * t125 + 0.1e1) * t92 + (t101 - t125) * qJD(3)) * t115) * t182) * t104) * t95 + (t90 * t170 + ((-t88 + t176) * t115 + (t152 * t101 + t153) * t116) * t95 * t136) * t104 * t126) * t134, 0, 0, 0; (-t108 * t111 + t112 * t189) * t169 + (-0.2e1 * t112 * t165 + (t111 * t94 - t112 * t93 + (t150 * t181 - (t151 * t133 + t156) * t125) * t113) * t109 + (t198 * t125 + t126 * t146) * t108) * t99, 0, -t160 * t200 + (t172 * t200 + (-0.2e1 * t145 * t201 + ((qJD(4) * t108 + 0.2e1 * t165) * t135 + (t135 * t93 + (-qJD(4) * t113 + t94) * t133) * t109) * t99) * t134) * t126 (t108 * t114 + t191) * t169 + (-0.2e1 * t166 + (-t109 * t114 + t108 - 0.2e1 * t190) * t94) * t99, 0, 0;];
JaD_rot  = t1;
