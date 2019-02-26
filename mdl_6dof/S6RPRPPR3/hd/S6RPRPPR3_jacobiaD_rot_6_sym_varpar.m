% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:24
% EndTime: 2019-02-26 20:40:25
% DurationCPUTime: 0.71s
% Computational Cost: add. (1787->91), mult. (2519->201), div. (480->12), fcn. (2968->9), ass. (0->94)
t125 = qJ(1) + pkin(9);
t123 = sin(t125);
t174 = qJD(3) * t123;
t121 = t123 ^ 2;
t132 = sin(qJ(3));
t127 = 0.1e1 / t132 ^ 2;
t134 = cos(qJ(3));
t130 = t134 ^ 2;
t180 = t127 * t130;
t117 = t121 * t180 + 0.1e1;
t115 = 0.1e1 / t117;
t126 = 0.1e1 / t132;
t161 = t127 * t174;
t124 = cos(t125);
t175 = qJD(1) * t134;
t162 = t124 * t175;
t173 = qJD(3) * t132;
t89 = ((t123 * t173 - t162) * t126 + t130 * t161) * t115;
t151 = -t89 + t174;
t152 = -t123 * t89 + qJD(3);
t131 = sin(qJ(6));
t133 = cos(qJ(6));
t150 = qJD(6) * t132 + qJD(1);
t172 = qJD(3) * t134;
t197 = t150 * t131 - t133 * t172;
t183 = t123 * t134;
t114 = atan2(-t183, t132);
t113 = cos(t114);
t112 = sin(t114);
t166 = t112 * t183;
t102 = t113 * t132 - t166;
t96 = 0.1e1 / t102;
t178 = t132 * t133;
t164 = t124 * t178;
t184 = t123 * t131;
t111 = t164 - t184;
t105 = 0.1e1 / t111;
t97 = 0.1e1 / t102 ^ 2;
t106 = 0.1e1 / t111 ^ 2;
t196 = t115 - 0.1e1;
t185 = t113 * t134;
t84 = t152 * t185 + (t151 * t132 - t162) * t112;
t195 = t84 * t96 * t97;
t122 = t124 ^ 2;
t94 = t122 * t130 * t97 + 0.1e1;
t92 = 0.1e1 / t94;
t194 = t92 * t97;
t193 = t96 * t92;
t179 = t131 * t132;
t110 = t123 * t133 + t124 * t179;
t104 = t110 ^ 2;
t103 = t104 * t106 + 0.1e1;
t187 = t106 * t110;
t149 = qJD(1) * t132 + qJD(6);
t145 = t149 * t133;
t91 = -t123 * t145 - t197 * t124;
t191 = t105 * t106 * t91;
t160 = t131 * t172;
t176 = qJD(1) * t124;
t90 = -qJD(6) * t164 - t124 * t160 - t133 * t176 + t149 * t184;
t192 = 0.1e1 / t103 ^ 2 * (-t104 * t191 - t90 * t187);
t189 = t124 * t97;
t186 = t112 * t132;
t182 = t124 * t131;
t181 = t126 * t130;
t177 = qJD(1) * t123;
t147 = t123 * t130 * t176;
t167 = t97 * t173;
t171 = 0.2e1 * (-t97 * t147 + (-t130 * t195 - t134 * t167) * t122) / t94 ^ 2;
t170 = 0.2e1 * t195;
t169 = 0.2e1 * t192;
t129 = t134 * t130;
t143 = qJD(3) * (-t127 * t129 - t134) * t126;
t168 = 0.2e1 / t117 ^ 2 * (t121 * t143 + t127 * t147);
t165 = t115 * t181;
t163 = t123 * t175;
t158 = t96 * t171;
t157 = t97 * t171;
t156 = 0.1e1 + t180;
t155 = 0.2e1 * t110 * t191;
t154 = t123 * t168;
t153 = t134 * t168;
t148 = t123 * t165;
t146 = t156 * t124;
t144 = -t105 * t131 + t133 * t187;
t99 = 0.1e1 / t103;
t142 = t144 * t99;
t109 = -t123 * t178 - t182;
t108 = -t123 * t179 + t124 * t133;
t101 = t156 * t123 * t115;
t88 = (t196 * t134 * t112 + t113 * t148) * t124;
t87 = t123 * t186 + t185 + (-t113 * t183 - t186) * t101;
t85 = -t156 * t154 + (qJD(1) * t146 + 0.2e1 * t123 * t143) * t115;
t1 = [t124 * t126 * t153 + (qJD(3) * t146 + t126 * t163) * t115, 0, t85, 0, 0, 0; (t173 * t193 + (t158 + (qJD(1) * t88 + t84) * t194) * t134) * t123 + (t88 * t157 * t134 + (t88 * t167 + (t88 * t170 + ((t89 * t148 + t196 * t173 + t153) * t112 + (t154 * t181 + t134 * t89 + (t129 * t161 - (t89 - 0.2e1 * t174) * t134) * t115) * t113) * t189) * t134 + (-t96 + (-(-t121 + t122) * t113 * t165 + t196 * t166) * t97) * t175) * t92) * t124, 0 (t177 * t193 + (t158 + (qJD(3) * t87 + t84) * t194) * t124) * t132 + (t87 * t124 * t157 + (-t124 * qJD(3) * t96 + (t124 * t170 + t97 * t177) * t87 + (-((-t101 * t176 - t123 * t85) * t113 + (-t152 * t101 + t151) * t112) * t134 - ((-t85 + t176) * t112 + (t151 * t101 - t152) * t113) * t132) * t189) * t92) * t134, 0, 0, 0; (-t105 * t108 + t109 * t187) * t169 + (t109 * t155 + (-t108 * t91 + t109 * t90 + (-t197 * t123 + t124 * t145) * t110) * t106 + (-t149 * t182 + (-t150 * t133 - t160) * t123) * t105) * t99, 0, t142 * t163 + (t142 * t173 + (t144 * t169 + ((qJD(6) * t105 + t155) * t133 + (t133 * t90 + (qJD(6) * t110 - t91) * t131) * t106) * t99) * t134) * t124, 0, 0, -0.2e1 * t192 + 0.2e1 * (-t90 * t106 * t99 + (-t106 * t192 - t99 * t191) * t110) * t110;];
JaD_rot  = t1;
