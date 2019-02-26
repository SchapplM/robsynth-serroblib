% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:12
% EndTime: 2019-02-26 21:03:13
% DurationCPUTime: 0.62s
% Computational Cost: add. (4795->72), mult. (2858->158), div. (686->14), fcn. (3330->7), ass. (0->75)
t123 = sin(qJ(1));
t117 = t123 ^ 2;
t115 = pkin(10) + qJ(3) + qJ(4);
t113 = sin(t115);
t108 = t113 ^ 2;
t114 = cos(t115);
t111 = 0.1e1 / t114 ^ 2;
t157 = t108 * t111;
t103 = t117 * t157 + 0.1e1;
t107 = t113 * t108;
t109 = t114 ^ 2;
t110 = 0.1e1 / t114;
t116 = qJD(3) + qJD(4);
t156 = t110 * t113;
t132 = t116 * (t107 * t110 / t109 + t156);
t124 = cos(qJ(1));
t148 = qJD(1) * t124;
t140 = t123 * t148;
t165 = 0.1e1 / t103 ^ 2 * (t117 * t132 + t140 * t157);
t173 = -0.2e1 * t165;
t101 = 0.1e1 / t103;
t138 = 0.1e1 + t157;
t171 = t123 * t138;
t96 = t101 * t171;
t172 = t123 * t96 - 0.1e1;
t152 = 0.1e1 / t123 * t124;
t122 = t124 ^ 2;
t170 = qJD(1) * (0.1e1 / t117 * t122 + 0.1e1) * t152;
t150 = t123 * t113;
t100 = atan2(-t150, -t114);
t98 = sin(t100);
t144 = t98 * t150;
t99 = cos(t100);
t95 = -t114 * t99 - t144;
t92 = 0.1e1 / t95;
t93 = 0.1e1 / t95 ^ 2;
t154 = t114 * t116;
t136 = t113 * t122 * t154;
t153 = t116 * t123;
t162 = t114 * t98;
t141 = t111 * t153;
t87 = (-(-t113 * t148 - t114 * t153) * t110 + t108 * t141) * t101;
t82 = (t87 - t153) * t162 + (-t98 * t148 + (-t123 * t87 + t116) * t99) * t113;
t168 = t82 * t92 * t93;
t90 = t108 * t122 * t93 + 0.1e1;
t169 = (t93 * t136 + (-t122 * t168 - t140 * t93) * t108) / t90 ^ 2;
t88 = 0.1e1 / t90;
t166 = t88 * t93;
t164 = t113 * t98;
t163 = t113 * t99;
t161 = t116 * t96;
t159 = t124 * t93;
t158 = t108 * t110;
t155 = t113 * t124;
t119 = 0.1e1 / t123 ^ 2;
t151 = t119 * t122;
t149 = qJD(1) * t123;
t147 = 0.2e1 * t168;
t106 = t109 * t151 + 0.1e1;
t146 = 0.2e1 / t106 ^ 2 * (-t109 * t170 - t119 * t136);
t145 = t92 * t169;
t143 = t88 * t154;
t142 = t123 * t158;
t139 = 0.2e1 * t93 * t169;
t137 = 0.1e1 + t151;
t135 = t138 * t124;
t134 = t137 * t113;
t131 = -t142 * t99 + t164;
t104 = 0.1e1 / t106;
t86 = (t101 * t131 - t164) * t124;
t85 = t113 * t146 * t152 + (qJD(1) * t134 - t152 * t154) * t104;
t84 = (-t123 + t96) * t162 - t172 * t163;
t83 = t171 * t173 + (qJD(1) * t135 + 0.2e1 * t123 * t132) * t101;
t80 = (-t92 * t88 * t149 + (-0.2e1 * t145 + (-t116 * t84 - t82) * t166) * t124) * t114 + (t84 * t124 * t139 + (-t124 * t116 * t92 - ((-t123 * t83 - t148 * t96) * t99 + (t172 * t87 + t153 - t161) * t98) * t93 * t155 + (t124 * t147 + t149 * t93) * t84 - ((t83 - t148) * t98 + (t87 * t96 + t116 + (-t87 - t161) * t123) * t99) * t114 * t159) * t88) * t113;
t1 = [t110 * t155 * t173 + (t116 * t135 - t149 * t156) * t101, 0, t83, t83, 0, 0; (-t92 * t143 + (0.2e1 * t145 + (qJD(1) * t86 + t82) * t166) * t113) * t123 + (-t86 * t93 * t143 + (t86 * t139 + (t86 * t147 + (t87 * t163 + t98 * t154 + 0.2e1 * t131 * t165 + ((-t142 * t87 - t154) * t98 + (t107 * t141 - (t87 - 0.2e1 * t153) * t113) * t99) * t101) * t159) * t88) * t113 + (-t92 + (-t144 + (t144 - (t117 - t122) * t99 * t158) * t101) * t93) * t113 * t88 * qJD(1)) * t124, 0, t80, t80, 0, 0; t137 * t114 * t146 + (0.2e1 * t114 * t170 + t116 * t134) * t104, 0, t85, t85, 0, 0;];
JaD_rot  = t1;
