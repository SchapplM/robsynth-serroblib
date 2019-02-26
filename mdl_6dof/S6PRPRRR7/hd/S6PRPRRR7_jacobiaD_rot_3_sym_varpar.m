% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR7_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiaD_rot_3_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:29
% EndTime: 2019-02-26 19:57:30
% DurationCPUTime: 0.43s
% Computational Cost: add. (1104->57), mult. (3515->153), div. (257->12), fcn. (4585->13), ass. (0->75)
t132 = sin(pkin(7));
t129 = t132 ^ 2;
t170 = 0.2e1 * t129;
t135 = cos(pkin(13));
t131 = sin(pkin(13));
t138 = sin(qJ(2));
t137 = cos(pkin(6));
t139 = cos(qJ(2));
t155 = t137 * t139;
t149 = -t131 * t138 + t135 * t155;
t133 = sin(pkin(6));
t136 = cos(pkin(7));
t159 = t133 * t136;
t115 = t149 * t132 + t135 * t159;
t160 = t132 * t139;
t124 = -t133 * t160 + t137 * t136;
t110 = atan2(t115, t124);
t103 = sin(t110);
t104 = cos(t110);
t92 = t103 * t115 + t104 * t124;
t89 = 0.1e1 / t92;
t130 = sin(pkin(14));
t134 = cos(pkin(14));
t156 = t137 * t138;
t147 = t131 * t156 - t135 * t139;
t148 = t131 * t155 + t135 * t138;
t150 = t131 * t132 * t133 - t136 * t148;
t102 = t150 * t130 - t134 * t147;
t98 = 0.1e1 / t102;
t121 = 0.1e1 / t124;
t99 = 0.1e1 / t102 ^ 2;
t122 = 0.1e1 / t124 ^ 2;
t90 = 0.1e1 / t92 ^ 2;
t116 = t131 * t159 + t132 * t148;
t114 = t116 ^ 2;
t120 = t147 * qJD(2);
t125 = -t131 * t139 - t135 * t156;
t118 = t125 * qJD(2);
t154 = qJD(2) * t133;
t152 = t138 * t154;
t158 = t133 * t138;
t162 = t115 * t122;
t151 = t158 * t162;
t113 = t115 ^ 2;
t109 = t113 * t122 + 0.1e1;
t105 = 0.1e1 / t109;
t163 = t105 * t132;
t84 = (-qJD(2) * t151 + t118 * t121) * t163;
t81 = (t115 * t84 + t132 * t152) * t104 + (t118 * t132 - t124 * t84) * t103;
t168 = t81 * t89 * t90;
t88 = t114 * t90 + 0.1e1;
t169 = (-t116 * t120 * t132 * t90 - t114 * t168) / t88 ^ 2;
t146 = -t121 * t125 + t151;
t85 = t146 * t163;
t167 = t115 * t85;
t166 = t124 * t85;
t119 = t148 * qJD(2);
t161 = t130 * t136;
t108 = -t119 * t134 + t120 * t161;
t165 = t98 * t99 * t108;
t101 = -t130 * t147 - t150 * t134;
t112 = -t134 * t148 + t147 * t161;
t164 = t101 * t112;
t157 = t134 * t136;
t153 = t122 * t129 * t138;
t123 = t121 * t122;
t117 = t149 * qJD(2);
t111 = -t130 * t148 - t147 * t157;
t107 = -t119 * t130 - t120 * t157;
t97 = t101 ^ 2;
t96 = t97 * t99 + 0.1e1;
t86 = 0.1e1 / t88;
t82 = (t132 * t158 - t167) * t104 + (t125 * t132 + t166) * t103;
t80 = t146 * (-t113 * t123 * t152 + t118 * t162) / t109 ^ 2 * t170 + (-t117 * t121 * t132 + (-t118 * t153 + (-t125 * t153 + (t123 * t133 * t138 ^ 2 * t170 - t122 * t160) * t115) * qJD(2)) * t133) * t105;
t1 = [0, t80, 0, 0, 0, 0; 0 (-((t115 * t80 + t84 * t166) * t104 + (-t124 * t80 + t84 * t167) * t103) * t90 * t86 + 0.2e1 * (t86 * t168 + t90 * t169) * t82) * t116 + (0.2e1 * t147 * t89 * t169 + (-t119 * t89 + (t82 * t120 + t147 * t81 + (-(-t118 * t85 + t125 * t84 + t139 * t154) * t104 - (-t117 + (qJD(2) * t85 - t84) * t158) * t103) * t116) * t90) * t86) * t132, 0, 0, 0, 0; 0, 0.2e1 * (-t111 * t98 + t99 * t164) / t96 ^ 2 * (t101 * t107 * t99 - t97 * t165) + ((-t119 * t157 + t120 * t130) * t98 + 0.2e1 * t164 * t165 + (-t111 * t108 - (t119 * t161 + t120 * t134) * t101 - t112 * t107) * t99) / t96, 0, 0, 0, 0;];
JaD_rot  = t1;
