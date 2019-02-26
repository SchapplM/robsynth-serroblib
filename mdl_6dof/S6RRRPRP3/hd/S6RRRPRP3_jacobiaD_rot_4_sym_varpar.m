% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:35
% EndTime: 2019-02-26 22:10:36
% DurationCPUTime: 0.70s
% Computational Cost: add. (3413->87), mult. (3414->191), div. (723->12), fcn. (4027->9), ass. (0->91)
t151 = sin(qJ(1));
t206 = 0.2e1 * t151;
t148 = qJ(2) + qJ(3);
t144 = cos(t148);
t150 = cos(pkin(10));
t180 = t151 * t150;
t149 = sin(pkin(10));
t152 = cos(qJ(1));
t184 = t149 * t152;
t130 = t144 * t184 - t180;
t124 = t130 ^ 2;
t181 = t151 * t149;
t183 = t150 * t152;
t131 = t144 * t183 + t181;
t126 = 0.1e1 / t131 ^ 2;
t121 = t124 * t126 + 0.1e1;
t128 = -t144 * t181 - t183;
t143 = sin(t148);
t145 = qJD(2) + qJD(3);
t185 = t145 * t152;
t170 = t143 * t185;
t122 = t128 * qJD(1) - t149 * t170;
t196 = t126 * t130;
t129 = -t144 * t180 + t184;
t123 = t129 * qJD(1) - t150 * t170;
t125 = 0.1e1 / t131;
t198 = t123 * t125 * t126;
t205 = (t122 * t196 - t124 * t198) / t121 ^ 2;
t146 = t151 ^ 2;
t139 = t143 ^ 2;
t141 = 0.1e1 / t144 ^ 2;
t192 = t139 * t141;
t137 = t146 * t192 + 0.1e1;
t135 = 0.1e1 / t137;
t140 = 0.1e1 / t144;
t178 = qJD(1) * t152;
t169 = t143 * t178;
t186 = t145 * t151;
t172 = t141 * t186;
t109 = (-(-t144 * t186 - t169) * t140 + t139 * t172) * t135;
t204 = t109 - t186;
t182 = t151 * t143;
t134 = atan2(-t182, -t144);
t133 = cos(t134);
t132 = sin(t134);
t173 = t132 * t182;
t117 = -t133 * t144 - t173;
t114 = 0.1e1 / t117;
t115 = 0.1e1 / t117 ^ 2;
t203 = t135 - 0.1e1;
t194 = t133 * t143;
t104 = (-t109 * t151 + t145) * t194 + (t144 * t204 - t169) * t132;
t202 = t104 * t114 * t115;
t138 = t143 * t139;
t189 = t140 * t143;
t160 = t145 * (t138 * t140 * t141 + t189);
t190 = t139 * t151;
t163 = t178 * t190;
t201 = (t141 * t163 + t146 * t160) / t137 ^ 2;
t200 = t115 * t143;
t199 = t115 * t152;
t197 = t125 * t149;
t195 = t132 * t151;
t193 = t139 * t140;
t147 = t152 ^ 2;
t191 = t139 * t147;
t188 = t143 * t152;
t187 = t144 * t145;
t179 = qJD(1) * t151;
t112 = t115 * t191 + 0.1e1;
t177 = 0.2e1 * (-t191 * t202 + (t143 * t147 * t187 - t163) * t115) / t112 ^ 2;
t176 = 0.2e1 * t202;
t175 = t115 * t188;
t174 = t130 * t198;
t171 = t145 * t182;
t168 = 0.1e1 + t192;
t167 = t143 * t177;
t166 = -0.2e1 * t143 * t201;
t165 = t201 * t206;
t164 = t133 * t135 * t193;
t162 = t168 * t152;
t161 = t150 * t196 - t197;
t119 = 0.1e1 / t121;
t118 = t168 * t151 * t135;
t110 = 0.1e1 / t112;
t108 = (t203 * t143 * t132 - t151 * t164) * t152;
t106 = -t144 * t195 + t194 + (t132 * t144 - t133 * t182) * t118;
t105 = -t168 * t165 + (qJD(1) * t162 + t160 * t206) * t135;
t102 = -0.2e1 * t161 * t188 * t205 + (t161 * t144 * t185 + (-0.2e1 * t174 * t183 + t179 * t197 + (t123 * t184 + (t122 * t152 - t130 * t179) * t150) * t126) * t143) * t119;
t101 = (t106 * t200 - t114 * t144) * t152 * t177 + ((-t114 * t179 + (-t106 * t145 - t104) * t199) * t144 + (-t114 * t185 - (-t105 * t133 * t151 - t204 * t132 + (t109 * t195 - t132 * t145 - t133 * t178) * t118) * t175 + (t115 * t179 + t152 * t176) * t106 - ((t105 - t178) * t132 + ((-t118 * t151 + 0.1e1) * t145 + (t118 - t151) * t109) * t133) * t144 * t199) * t143) * t110;
t1 = [t140 * t152 * t166 + (t145 * t162 - t179 * t189) * t135, t105, t105, 0, 0, 0; (t114 * t167 + (-t114 * t187 + (qJD(1) * t108 + t104) * t200) * t110) * t151 + (t115 * t167 * t108 + (-((t166 - t187 + (t109 * t140 * t190 + t187) * t135) * t132 + (t165 * t193 - t109 * t143 + (-t138 * t172 + (t109 - 0.2e1 * t186) * t143) * t135) * t133) * t175 + (-t115 * t187 + t143 * t176) * t108 + (-t114 + ((-t146 + t147) * t164 + t203 * t173) * t115) * t143 * qJD(1)) * t110) * t152, t101, t101, 0, 0, 0; 0.2e1 * (-t125 * t128 + t129 * t196) * t205 + ((-t130 * qJD(1) + t149 * t171) * t125 + 0.2e1 * t129 * t174 + (-t128 * t123 - (-t131 * qJD(1) + t150 * t171) * t130 - t129 * t122) * t126) * t119, t102, t102, 0, 0, 0;];
JaD_rot  = t1;
