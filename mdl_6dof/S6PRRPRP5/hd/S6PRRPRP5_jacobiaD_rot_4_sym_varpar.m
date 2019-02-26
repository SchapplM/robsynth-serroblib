% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRP5_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:31
% EndTime: 2019-02-26 20:03:32
% DurationCPUTime: 0.72s
% Computational Cost: add. (2404->89), mult. (7372->192), div. (524->12), fcn. (9540->11), ass. (0->86)
t161 = sin(pkin(10));
t163 = cos(pkin(10));
t166 = sin(qJ(2));
t164 = cos(pkin(6));
t168 = cos(qJ(2));
t189 = t164 * t168;
t151 = -t161 * t166 + t163 * t189;
t141 = t151 * qJD(2);
t190 = t164 * t166;
t152 = t161 * t168 + t163 * t190;
t165 = sin(qJ(3));
t162 = sin(pkin(6));
t193 = t162 * t165;
t181 = t163 * t193;
t167 = cos(qJ(3));
t186 = qJD(3) * t167;
t118 = -qJD(3) * t181 + t141 * t165 + t152 * t186;
t192 = t162 * t167;
t134 = t152 * t165 + t163 * t192;
t131 = t134 ^ 2;
t155 = -t164 * t167 + t166 * t193;
t149 = 0.1e1 / t155 ^ 2;
t129 = t131 * t149 + 0.1e1;
t126 = 0.1e1 / t129;
t156 = t164 * t165 + t166 * t192;
t187 = qJD(2) * t168;
t180 = t162 * t187;
t139 = t156 * qJD(3) + t165 * t180;
t148 = 0.1e1 / t155;
t197 = t134 * t149;
t106 = (-t118 * t148 + t139 * t197) * t126;
t130 = atan2(-t134, t155);
t122 = sin(t130);
t123 = cos(t130);
t178 = -t122 * t155 - t123 * t134;
t103 = t178 * t106 - t122 * t118 + t123 * t139;
t117 = -t122 * t134 + t123 * t155;
t114 = 0.1e1 / t117;
t115 = 0.1e1 / t117 ^ 2;
t205 = t103 * t114 * t115;
t182 = t161 * t190;
t154 = t163 * t168 - t182;
t176 = -t154 * t165 + t161 * t192;
t204 = -0.2e1 * t176 * t205;
t191 = t162 * t168;
t175 = -t148 * t151 + t191 * t197;
t203 = t165 * t175;
t195 = t139 * t148 * t149;
t202 = -0.2e1 * (t118 * t197 - t131 * t195) / t129 ^ 2;
t153 = t161 * t189 + t163 * t166;
t145 = 0.1e1 / t153;
t146 = 0.1e1 / t153 ^ 2;
t201 = t115 * t176;
t138 = t154 * t167 + t161 * t193;
t143 = t153 * qJD(2);
t120 = t138 * qJD(3) - t143 * t165;
t200 = t120 * t115;
t199 = t122 * t176;
t198 = t123 * t176;
t196 = t138 * t154;
t194 = t153 * t165;
t188 = qJD(2) * t166;
t132 = t176 ^ 2;
t112 = t132 * t115 + 0.1e1;
t185 = 0.2e1 * (-t132 * t205 - t176 * t200) / t112 ^ 2;
t121 = t176 * qJD(3) - t143 * t167;
t133 = t138 ^ 2;
t128 = t133 * t146 + 0.1e1;
t144 = -qJD(2) * t182 + t163 * t187;
t147 = t145 * t146;
t184 = 0.2e1 * (t138 * t146 * t121 - t133 * t147 * t144) / t128 ^ 2;
t179 = -0.2e1 * t134 * t195;
t136 = t152 * t167 - t181;
t177 = -t136 * t148 + t156 * t197;
t142 = t152 * qJD(2);
t140 = -t155 * qJD(3) + t167 * t180;
t124 = 0.1e1 / t128;
t119 = -t134 * qJD(3) + t141 * t167;
t109 = 0.1e1 / t112;
t108 = t126 * t203;
t107 = t177 * t126;
t105 = (-t122 * t151 + t123 * t191) * t165 + t178 * t108;
t104 = t178 * t107 - t122 * t136 + t123 * t156;
t102 = t177 * t202 + (t156 * t179 - t119 * t148 + (t118 * t156 + t134 * t140 + t136 * t139) * t149) * t126;
t100 = t202 * t203 + (t175 * t186 + (t179 * t191 + t142 * t148 + (t139 * t151 + (t118 * t168 - t134 * t188) * t162) * t149) * t165) * t126;
t1 = [0, t100, t102, 0, 0, 0; 0 (-t105 * t201 + t114 * t194) * t185 + ((-t144 * t165 - t153 * t186) * t114 + (-t200 + t204) * t105 + (t194 * t103 + (-t100 * t134 - t108 * t118 + (-t165 * t188 + t168 * t186) * t162 + (-t108 * t155 - t151 * t165) * t106) * t198 + (-t151 * t186 - t100 * t155 - t108 * t139 + t142 * t165 + (t108 * t134 - t165 * t191) * t106) * t199) * t115) * t109 (-t104 * t201 - t114 * t138) * t185 + (t104 * t204 + t121 * t114 + (-t138 * t103 - t104 * t120 + (-t102 * t134 - t107 * t118 + t140 + (-t107 * t155 - t136) * t106) * t198 + (-t102 * t155 - t107 * t139 - t119 + (t107 * t134 - t156) * t106) * t199) * t115) * t109, 0, 0, 0; 0 (t145 * t153 * t167 + t146 * t196) * t184 + (qJD(3) * t145 * t194 + (-t121 * t154 + t138 * t143) * t146 + (0.2e1 * t147 * t196 + (t146 * t153 - t145) * t167) * t144) * t124, -t176 * t145 * t184 + (-t144 * t146 * t176 - t120 * t145) * t124, 0, 0, 0;];
JaD_rot  = t1;
