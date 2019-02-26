% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:47
% EndTime: 2019-02-26 19:49:48
% DurationCPUTime: 0.68s
% Computational Cost: add. (2385->85), mult. (7313->189), div. (516->12), fcn. (9465->11), ass. (0->86)
t158 = sin(pkin(10));
t160 = cos(pkin(10));
t165 = cos(qJ(2));
t161 = cos(pkin(6));
t163 = sin(qJ(2));
t185 = t161 * t163;
t153 = t158 * t165 + t160 * t185;
t144 = t153 * qJD(2);
t159 = sin(pkin(6));
t164 = cos(qJ(4));
t162 = sin(qJ(4));
t184 = t161 * t165;
t174 = -t158 * t163 + t160 * t184;
t172 = t174 * t162;
t121 = qJD(4) * t172 + (t160 * t159 * qJD(4) + t144) * t164;
t189 = t159 * t162;
t136 = t160 * t189 - t174 * t164;
t133 = t136 ^ 2;
t186 = t159 * t165;
t156 = t161 * t162 + t164 * t186;
t151 = 0.1e1 / t156 ^ 2;
t129 = t133 * t151 + 0.1e1;
t126 = 0.1e1 / t129;
t157 = t161 * t164 - t162 * t186;
t188 = t159 * t163;
t177 = qJD(2) * t188;
t139 = t157 * qJD(4) - t164 * t177;
t150 = 0.1e1 / t156;
t195 = t136 * t151;
t106 = (t121 * t150 - t139 * t195) * t126;
t130 = atan2(t136, t156);
t122 = sin(t130);
t123 = cos(t130);
t176 = -t122 * t156 + t123 * t136;
t103 = t176 * t106 + t122 * t121 + t123 * t139;
t117 = t122 * t136 + t123 * t156;
t114 = 0.1e1 / t117;
t115 = 0.1e1 / t117 ^ 2;
t203 = t103 * t114 * t115;
t154 = t158 * t184 + t160 * t163;
t134 = -t154 * t164 + t158 * t189;
t202 = 0.2e1 * t134 * t203;
t193 = t139 * t150 * t151;
t201 = (t121 * t195 - t133 * t193) / t129 ^ 2;
t145 = t154 * qJD(2);
t155 = -t158 * t185 + t160 * t165;
t148 = 0.1e1 / t155 ^ 2;
t200 = t145 * t148;
t178 = t136 * t188;
t173 = t150 * t153 + t151 * t178;
t199 = t164 * t173;
t147 = 0.1e1 / t155;
t198 = t115 * t134;
t197 = t122 * t134;
t196 = t123 * t134;
t194 = t136 * t157;
t192 = t147 * t200;
t191 = t148 * t154;
t190 = t155 * t164;
t187 = t159 * t164;
t183 = qJD(2) * t165;
t182 = qJD(4) * t162;
t131 = t134 ^ 2;
t112 = t131 * t115 + 0.1e1;
t135 = t154 * t162 + t158 * t187;
t146 = t155 * qJD(2);
t119 = t135 * qJD(4) - t146 * t164;
t181 = 0.2e1 * (t119 * t198 - t131 * t203) / t112 ^ 2;
t118 = -t134 * qJD(4) + t146 * t162;
t132 = t135 ^ 2;
t128 = t132 * t148 + 0.1e1;
t180 = 0.2e1 * (t135 * t148 * t118 + t132 * t192) / t128 ^ 2;
t137 = t160 * t187 + t172;
t175 = -t137 * t150 + t151 * t194;
t143 = t174 * qJD(2);
t138 = -t156 * qJD(4) + t162 * t177;
t124 = 0.1e1 / t128;
t120 = t136 * qJD(4) + t144 * t162;
t109 = 0.1e1 / t112;
t108 = t126 * t199;
t107 = t175 * t126;
t105 = (t122 * t153 - t123 * t188) * t164 + t176 * t108;
t104 = -t176 * t107 + t122 * t137 + t123 * t157;
t102 = 0.2e1 * t175 * t201 + (0.2e1 * t193 * t194 - t120 * t150 + (-t121 * t157 - t136 * t138 - t137 * t139) * t151) * t126;
t100 = -0.2e1 * t199 * t201 + (-t173 * t182 + (-0.2e1 * t178 * t193 + t143 * t150 + (-t139 * t153 + (t121 * t163 + t136 * t183) * t159) * t151) * t164) * t126;
t1 = [0, t100, 0, t102, 0, 0; 0 (t105 * t198 + t114 * t190) * t181 + ((t145 * t164 + t155 * t182) * t114 + t105 * t202 + (-t105 * t119 + t190 * t103 - (t100 * t136 + t108 * t121 + (t163 * t182 - t164 * t183) * t159 + (-t108 * t156 + t153 * t164) * t106) * t196 - (-t153 * t182 - t100 * t156 - t108 * t139 + t143 * t164 + (-t108 * t136 + t163 * t187) * t106) * t197) * t115) * t109, 0 (t104 * t198 - t114 * t135) * t181 + (t104 * t202 + t118 * t114 + (-t135 * t103 - t104 * t119 - (t102 * t136 - t107 * t121 + t138 + (t107 * t156 + t137) * t106) * t196 - (-t102 * t156 + t107 * t139 - t120 + (t107 * t136 - t157) * t106) * t197) * t115) * t109, 0, 0; 0 (-t135 * t191 - t162) * t180 + (t118 * t191 + qJD(4) * t164 + (t146 * t148 + 0.2e1 * t154 * t192) * t135) * t124, 0, t134 * t147 * t180 + (-t119 * t147 - t134 * t200) * t124, 0, 0;];
JaD_rot  = t1;
