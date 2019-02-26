% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP4_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:01
% EndTime: 2019-02-26 19:52:02
% DurationCPUTime: 0.40s
% Computational Cost: add. (941->56), mult. (2271->131), div. (423->14), fcn. (2956->11), ass. (0->66)
t146 = sin(qJ(2));
t147 = cos(qJ(2));
t144 = sin(pkin(10));
t173 = cos(pkin(6));
t161 = t144 * t173;
t172 = cos(pkin(10));
t132 = -t146 * t161 + t172 * t147;
t140 = pkin(11) + qJ(4);
t136 = sin(t140);
t137 = cos(t140);
t145 = sin(pkin(6));
t166 = t144 * t145;
t155 = -t132 * t136 + t137 * t166;
t177 = t155 * qJD(4);
t158 = t173 * t172;
t128 = t144 * t146 - t147 * t158;
t165 = t145 * t147;
t118 = atan2(-t128, -t165);
t116 = sin(t118);
t117 = cos(t118);
t103 = -t116 * t128 - t117 * t165;
t100 = 0.1e1 / t103;
t115 = t132 * t137 + t136 * t166;
t111 = 0.1e1 / t115;
t141 = 0.1e1 / t147;
t101 = 0.1e1 / t103 ^ 2;
t112 = 0.1e1 / t115 ^ 2;
t142 = 0.1e1 / t147 ^ 2;
t130 = t144 * t147 + t146 * t158;
t123 = t130 * qJD(2);
t164 = qJD(2) * t146;
t167 = t142 * t146;
t162 = t128 * t167;
t126 = t128 ^ 2;
t139 = 0.1e1 / t145 ^ 2;
t121 = t126 * t139 * t142 + 0.1e1;
t119 = 0.1e1 / t121;
t138 = 0.1e1 / t145;
t168 = t119 * t138;
t95 = (qJD(2) * t162 + t123 * t141) * t168;
t92 = (-t128 * t95 + t145 * t164) * t117 + (t95 * t165 - t123) * t116;
t176 = t100 * t101 * t92;
t110 = t155 ^ 2;
t106 = t110 * t112 + 0.1e1;
t154 = -t172 * t146 - t147 * t161;
t124 = t154 * qJD(2);
t108 = t115 * qJD(4) + t124 * t136;
t169 = t112 * t155;
t109 = t124 * t137 + t177;
t170 = t109 * t111 * t112;
t175 = 0.1e1 / t106 ^ 2 * (-t108 * t169 - t110 * t170);
t156 = t130 * t141 + t162;
t96 = t156 * t168;
t174 = t128 * t96;
t171 = t101 * t154;
t163 = -0.2e1 * t175;
t157 = -t111 * t136 - t137 * t169;
t143 = t141 * t142;
t127 = t154 ^ 2;
t125 = t132 * qJD(2);
t122 = t128 * qJD(2);
t104 = 0.1e1 / t106;
t99 = t127 * t101 + 0.1e1;
t93 = (t145 * t146 - t174) * t117 + (t96 * t165 - t130) * t116;
t91 = (-0.2e1 * t156 / t121 ^ 2 * (t123 * t128 * t142 + t126 * t143 * t164) * t139 + (t123 * t167 - t122 * t141 + (t130 * t167 + (0.2e1 * t143 * t146 ^ 2 + t141) * t128) * qJD(2)) * t119) * t138;
t1 = [0, t91, 0, 0, 0, 0; 0, 0.2e1 * (-t100 * t132 - t93 * t171) / t99 ^ 2 * (-t125 * t171 - t127 * t176) + (t124 * t100 + (-t93 * t125 - t132 * t92) * t101 - (0.2e1 * t93 * t176 + (-(-t123 * t96 - t128 * t91 - t130 * t95 + (t95 * t96 + qJD(2)) * t165) * t117 - (t95 * t174 + t122 + (t147 * t91 + (-qJD(2) * t96 - t95) * t146) * t145) * t116) * t101) * t154) / t99, 0, 0, 0, 0; 0, -t157 * t154 * t163 + (t157 * t125 - ((-qJD(4) * t111 + 0.2e1 * t155 * t170) * t137 + (t108 * t137 + (t109 + t177) * t136) * t112) * t154) * t104, 0, t163 - 0.2e1 * (t104 * t108 * t112 - (-t104 * t170 - t112 * t175) * t155) * t155, 0, 0;];
JaD_rot  = t1;
