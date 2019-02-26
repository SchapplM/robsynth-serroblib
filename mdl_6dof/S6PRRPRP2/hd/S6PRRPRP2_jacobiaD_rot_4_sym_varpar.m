% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRP2_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:48
% EndTime: 2019-02-26 20:01:48
% DurationCPUTime: 0.38s
% Computational Cost: add. (941->56), mult. (2271->131), div. (423->14), fcn. (2956->11), ass. (0->66)
t149 = sin(qJ(2));
t150 = cos(qJ(2));
t147 = sin(pkin(10));
t176 = cos(pkin(6));
t164 = t147 * t176;
t175 = cos(pkin(10));
t135 = -t149 * t164 + t175 * t150;
t143 = qJ(3) + pkin(11);
t139 = sin(t143);
t140 = cos(t143);
t148 = sin(pkin(6));
t169 = t147 * t148;
t158 = -t135 * t139 + t140 * t169;
t180 = t158 * qJD(3);
t161 = t176 * t175;
t131 = t147 * t149 - t150 * t161;
t168 = t148 * t150;
t121 = atan2(-t131, -t168);
t119 = sin(t121);
t120 = cos(t121);
t106 = -t119 * t131 - t120 * t168;
t103 = 0.1e1 / t106;
t118 = t135 * t140 + t139 * t169;
t114 = 0.1e1 / t118;
t144 = 0.1e1 / t150;
t104 = 0.1e1 / t106 ^ 2;
t115 = 0.1e1 / t118 ^ 2;
t145 = 0.1e1 / t150 ^ 2;
t133 = t147 * t150 + t149 * t161;
t126 = t133 * qJD(2);
t167 = qJD(2) * t149;
t170 = t145 * t149;
t165 = t131 * t170;
t129 = t131 ^ 2;
t142 = 0.1e1 / t148 ^ 2;
t124 = t129 * t142 * t145 + 0.1e1;
t122 = 0.1e1 / t124;
t141 = 0.1e1 / t148;
t171 = t122 * t141;
t98 = (qJD(2) * t165 + t126 * t144) * t171;
t95 = (-t131 * t98 + t148 * t167) * t120 + (t98 * t168 - t126) * t119;
t179 = t103 * t104 * t95;
t113 = t158 ^ 2;
t109 = t113 * t115 + 0.1e1;
t157 = -t175 * t149 - t150 * t164;
t127 = t157 * qJD(2);
t111 = t118 * qJD(3) + t127 * t139;
t172 = t115 * t158;
t112 = t127 * t140 + t180;
t173 = t112 * t114 * t115;
t178 = 0.1e1 / t109 ^ 2 * (-t111 * t172 - t113 * t173);
t159 = t133 * t144 + t165;
t99 = t159 * t171;
t177 = t131 * t99;
t174 = t104 * t157;
t166 = -0.2e1 * t178;
t160 = -t114 * t139 - t140 * t172;
t146 = t144 * t145;
t130 = t157 ^ 2;
t128 = t135 * qJD(2);
t125 = t131 * qJD(2);
t107 = 0.1e1 / t109;
t102 = t104 * t130 + 0.1e1;
t96 = (t148 * t149 - t177) * t120 + (t99 * t168 - t133) * t119;
t94 = (-0.2e1 * t159 / t124 ^ 2 * (t126 * t131 * t145 + t129 * t146 * t167) * t142 + (t126 * t170 - t125 * t144 + (t133 * t170 + (0.2e1 * t146 * t149 ^ 2 + t144) * t131) * qJD(2)) * t122) * t141;
t1 = [0, t94, 0, 0, 0, 0; 0, 0.2e1 * (-t103 * t135 - t96 * t174) * (-t128 * t174 - t130 * t179) / t102 ^ 2 + (t127 * t103 + (-t96 * t128 - t135 * t95) * t104 - (0.2e1 * t96 * t179 + (-(-t126 * t99 - t131 * t94 - t133 * t98 + (t98 * t99 + qJD(2)) * t168) * t120 - (t98 * t177 + t125 + (t150 * t94 + (-qJD(2) * t99 - t98) * t149) * t148) * t119) * t104) * t157) / t102, 0, 0, 0, 0; 0, -t160 * t157 * t166 + (t160 * t128 - ((-qJD(3) * t114 + 0.2e1 * t158 * t173) * t140 + (t111 * t140 + (t112 + t180) * t139) * t115) * t157) * t107, t166 - 0.2e1 * (t107 * t111 * t115 - (-t107 * t173 - t115 * t178) * t158) * t158, 0, 0, 0;];
JaD_rot  = t1;
