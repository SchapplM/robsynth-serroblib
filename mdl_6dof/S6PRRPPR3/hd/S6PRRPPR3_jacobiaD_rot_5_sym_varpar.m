% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:23
% EndTime: 2019-02-26 19:59:24
% DurationCPUTime: 0.38s
% Computational Cost: add. (694->54), mult. (2307->132), div. (427->14), fcn. (2998->11), ass. (0->67)
t145 = sin(pkin(10));
t147 = cos(pkin(10));
t150 = sin(qJ(2));
t148 = cos(pkin(6));
t152 = cos(qJ(2));
t170 = t148 * t152;
t133 = t145 * t150 - t147 * t170;
t146 = sin(pkin(6));
t172 = t146 * t152;
t125 = atan2(t133, t172);
t121 = sin(t125);
t122 = cos(t125);
t108 = t121 * t133 + t122 * t172;
t105 = 0.1e1 / t108;
t149 = sin(qJ(3));
t151 = cos(qJ(3));
t171 = t148 * t150;
t159 = t145 * t171 - t147 * t152;
t173 = t145 * t146;
t161 = t149 * t159 + t151 * t173;
t115 = 0.1e1 / t161;
t142 = 0.1e1 / t152;
t106 = 0.1e1 / t108 ^ 2;
t116 = 0.1e1 / t161 ^ 2;
t143 = 0.1e1 / t152 ^ 2;
t134 = t145 * t152 + t147 * t171;
t128 = t134 * qJD(2);
t174 = t143 * t150;
t165 = t133 * t174;
t131 = t133 ^ 2;
t141 = 0.1e1 / t146 ^ 2;
t126 = t131 * t141 * t143 + 0.1e1;
t123 = 0.1e1 / t126;
t140 = 0.1e1 / t146;
t175 = t123 * t140;
t100 = (qJD(2) * t165 + t128 * t142) * t175;
t163 = -t121 * t172 + t122 * t133;
t166 = t122 * t146 * t150;
t97 = -qJD(2) * t166 + t163 * t100 + t121 * t128;
t180 = t105 * t106 * t97;
t160 = t145 * t170 + t147 * t150;
t179 = t106 * t160;
t120 = t149 * t173 - t151 * t159;
t118 = t120 ^ 2;
t178 = t116 * t118;
t177 = t116 * t120;
t117 = t115 * t116;
t176 = t117 * t118;
t169 = qJD(3) * t120;
t112 = 0.1e1 + t178;
t129 = t160 * qJD(2);
t113 = -t129 * t149 + t169;
t114 = t161 * qJD(3) - t129 * t151;
t167 = t114 * t177;
t168 = 0.2e1 / t112 ^ 2 * (t113 * t176 + t167);
t164 = t115 * t151 + t149 * t177;
t162 = t134 * t142 + t165;
t144 = t142 * t143;
t132 = t160 ^ 2;
t130 = t159 * qJD(2);
t127 = t133 * qJD(2);
t110 = 0.1e1 / t112;
t104 = t106 * t132 + 0.1e1;
t101 = t162 * t175;
t98 = t163 * t101 + t121 * t134 - t166;
t96 = (-0.2e1 * t162 / t126 ^ 2 * (qJD(2) * t131 * t144 * t150 + t128 * t133 * t143) * t141 + (t128 * t174 - t127 * t142 + (t134 * t174 + (0.2e1 * t144 * t150 ^ 2 + t142) * t133) * qJD(2)) * t123) * t140;
t1 = [0, t96, 0, 0, 0, 0; 0, 0.2e1 * (-t105 * t159 - t98 * t179) * (-t130 * t179 - t132 * t180) / t104 ^ 2 + (t129 * t105 + (-t98 * t130 - t159 * t97) * t106 - (0.2e1 * t98 * t180 + (-(-qJD(2) * t172 + t101 * t128 + t133 * t96 + (-t101 * t172 + t134) * t100) * t122 - (-t100 * t101 * t133 - t127 + (-t152 * t96 + (qJD(2) * t101 + t100) * t150) * t146) * t121) * t106) * t160) / t104, 0, 0, 0, 0; 0, -t164 * t160 * t168 + (-t164 * t130 - ((-0.2e1 * t113 * t117 * t120 + qJD(3) * t115) * t149 + (-t114 * t149 + (-t113 - t169) * t151) * t116) * t160) * t110 (t115 * t161 + t178) * t168 + (-0.2e1 * t167 + (-t116 * t161 + t115 - 0.2e1 * t176) * t113) * t110, 0, 0, 0;];
JaD_rot  = t1;
