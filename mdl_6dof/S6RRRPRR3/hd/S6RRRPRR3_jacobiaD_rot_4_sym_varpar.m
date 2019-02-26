% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR3_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:07
% EndTime: 2019-02-26 22:17:08
% DurationCPUTime: 0.61s
% Computational Cost: add. (3098->71), mult. (2860->161), div. (736->13), fcn. (3352->7), ass. (0->75)
t117 = qJ(2) + qJ(3);
t109 = cos(t117);
t105 = 0.1e1 / t109;
t169 = -0.2e1 * t105;
t118 = sin(qJ(1));
t108 = sin(t117);
t104 = t108 ^ 2;
t106 = 0.1e1 / t109 ^ 2;
t152 = t104 * t106;
t134 = 0.1e1 + t152;
t112 = t118 ^ 2;
t99 = t112 * t152 + 0.1e1;
t97 = 0.1e1 / t99;
t130 = t134 * t97;
t92 = t118 * t130;
t168 = t118 * t92 - 0.1e1;
t119 = cos(qJ(1));
t114 = 0.1e1 / t119;
t148 = t114 * t118;
t113 = t119 ^ 2;
t167 = qJD(1) * (t112 / t113 + 0.1e1) * t148;
t147 = t118 * t108;
t96 = atan2(-t147, -t109);
t94 = sin(t96);
t139 = t94 * t147;
t95 = cos(t96);
t91 = -t109 * t95 - t139;
t88 = 0.1e1 / t91;
t89 = 0.1e1 / t91 ^ 2;
t166 = 0.2e1 * t108;
t165 = t97 - 0.1e1;
t144 = qJD(1) * t119;
t136 = t118 * t144;
t110 = qJD(2) + qJD(3);
t151 = t109 * t110;
t138 = t89 * t151;
t150 = t110 * t118;
t156 = t109 * t94;
t83 = (-(-t108 * t144 - t109 * t150) * t105 + t150 * t152) * t97;
t78 = (t83 - t150) * t156 + (-t94 * t144 + (-t118 * t83 + t110) * t95) * t108;
t163 = t78 * t88 * t89;
t158 = t104 * t89;
t86 = t113 * t158 + 0.1e1;
t164 = (t113 * t108 * t138 + (-t113 * t163 - t136 * t89) * t104) / t86 ^ 2;
t84 = 0.1e1 / t86;
t161 = t84 * t89;
t103 = t108 * t104;
t107 = t105 * t106;
t127 = t110 * (t103 * t107 + t105 * t108);
t160 = (t112 * t127 + t136 * t152) / t99 ^ 2;
t159 = t88 * t84;
t157 = t105 * t97;
t155 = t110 * t92;
t153 = t119 * t89;
t149 = t112 / t119 ^ 2;
t146 = qJD(1) * t108;
t145 = qJD(1) * t118;
t143 = 0.2e1 * t163;
t102 = t106 * t149 + 0.1e1;
t142 = 0.2e1 / t102 ^ 2 * (t107 * t108 * t110 * t149 + t106 * t167);
t141 = t88 * t164;
t140 = t118 * t157;
t137 = t165 * t108;
t135 = 0.2e1 * t89 * t164;
t133 = 0.1e1 + t149;
t132 = t160 * t169;
t131 = t104 * t140;
t128 = t133 * t108 * t106;
t100 = 0.1e1 / t102;
t82 = (-t131 * t95 + t137 * t94) * t119;
t81 = t106 * t114 * t142 * t147 + ((-0.2e1 * t104 * t107 - t105) * t110 * t148 - qJD(1) * t128) * t100;
t80 = (-t118 + t92) * t156 - t168 * t95 * t108;
t79 = t130 * t144 + 0.2e1 * (t127 * t97 - t134 * t160) * t118;
t76 = (-t145 * t159 + (-0.2e1 * t141 + (-t110 * t80 - t78) * t161) * t119) * t109 + (t80 * t119 * t135 + (-t119 * t110 * t88 + (t119 * t143 + t145 * t89) * t80 + (-((-t118 * t79 - t144 * t92) * t95 + (t168 * t83 + t150 - t155) * t94) * t108 - ((t79 - t144) * t94 + (t83 * t92 + t110 + (-t83 - t155) * t118) * t95) * t109) * t153) * t84) * t108;
t1 = [-t140 * t146 + (t108 * t132 + t110 * t130) * t119, t79, t79, 0, 0, 0; (-t151 * t159 + (0.2e1 * t141 + (qJD(1) * t82 + t78) * t161) * t108) * t118 + (t82 * t135 * t108 + (-t82 * t138 + (t82 * t143 + ((-t131 * t83 - t151 * t165 + t160 * t166) * t94 + (-t83 * t137 + (t104 * t132 + (t103 * t106 + t166) * t97 * t110) * t118) * t95) * t153) * t108 + (-t88 + t165 * t89 * t139 - (t112 - t113) * t95 * t157 * t158) * t146) * t84) * t119, t76, t76, 0, 0, 0; t133 * t105 * t142 + (-t110 * t128 + t167 * t169) * t100, t81, t81, 0, 0, 0;];
JaD_rot  = t1;
