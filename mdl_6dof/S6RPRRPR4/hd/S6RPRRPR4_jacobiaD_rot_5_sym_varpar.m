% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:47
% EndTime: 2019-02-26 21:02:48
% DurationCPUTime: 0.66s
% Computational Cost: add. (5265->87), mult. (3414->191), div. (723->12), fcn. (4027->9), ass. (0->91)
t152 = sin(qJ(1));
t207 = 0.2e1 * t152;
t146 = pkin(10) + qJ(3) + qJ(4);
t145 = cos(t146);
t151 = cos(pkin(11));
t181 = t152 * t151;
t150 = sin(pkin(11));
t153 = cos(qJ(1));
t185 = t150 * t153;
t133 = t145 * t185 - t181;
t127 = t133 ^ 2;
t182 = t152 * t150;
t184 = t151 * t153;
t134 = t145 * t184 + t182;
t129 = 0.1e1 / t134 ^ 2;
t122 = t127 * t129 + 0.1e1;
t131 = -t145 * t182 - t184;
t144 = sin(t146);
t147 = qJD(3) + qJD(4);
t186 = t147 * t153;
t171 = t144 * t186;
t123 = t131 * qJD(1) - t150 * t171;
t195 = t129 * t133;
t132 = -t145 * t181 + t185;
t124 = t132 * qJD(1) - t151 * t171;
t128 = 0.1e1 / t134;
t199 = t124 * t128 * t129;
t206 = (t123 * t195 - t127 * t199) / t122 ^ 2;
t148 = t152 ^ 2;
t140 = t144 ^ 2;
t142 = 0.1e1 / t145 ^ 2;
t193 = t140 * t142;
t138 = t148 * t193 + 0.1e1;
t136 = 0.1e1 / t138;
t141 = 0.1e1 / t145;
t179 = qJD(1) * t153;
t170 = t144 * t179;
t187 = t147 * t152;
t173 = t142 * t187;
t110 = (-(-t145 * t187 - t170) * t141 + t140 * t173) * t136;
t205 = t110 - t187;
t183 = t152 * t144;
t135 = atan2(-t183, -t145);
t126 = cos(t135);
t125 = sin(t135);
t174 = t125 * t183;
t118 = -t126 * t145 - t174;
t115 = 0.1e1 / t118;
t116 = 0.1e1 / t118 ^ 2;
t204 = t136 - 0.1e1;
t197 = t126 * t144;
t105 = (-t110 * t152 + t147) * t197 + (t205 * t145 - t170) * t125;
t203 = t105 * t115 * t116;
t139 = t144 * t140;
t190 = t141 * t144;
t161 = t147 * (t139 * t141 * t142 + t190);
t191 = t140 * t152;
t164 = t179 * t191;
t202 = (t142 * t164 + t148 * t161) / t138 ^ 2;
t201 = t116 * t144;
t200 = t116 * t153;
t198 = t125 * t152;
t196 = t128 * t150;
t194 = t140 * t141;
t149 = t153 ^ 2;
t192 = t140 * t149;
t189 = t144 * t153;
t188 = t145 * t147;
t180 = qJD(1) * t152;
t113 = t116 * t192 + 0.1e1;
t178 = 0.2e1 * (-t192 * t203 + (t144 * t149 * t188 - t164) * t116) / t113 ^ 2;
t177 = 0.2e1 * t203;
t176 = t116 * t189;
t175 = t133 * t199;
t172 = t147 * t183;
t169 = 0.1e1 + t193;
t168 = t144 * t178;
t167 = -0.2e1 * t144 * t202;
t166 = t202 * t207;
t165 = t126 * t136 * t194;
t163 = t169 * t153;
t162 = t151 * t195 - t196;
t120 = 0.1e1 / t122;
t119 = t169 * t152 * t136;
t111 = 0.1e1 / t113;
t109 = (t204 * t144 * t125 - t152 * t165) * t153;
t107 = -t145 * t198 + t197 + (t125 * t145 - t126 * t183) * t119;
t106 = -t169 * t166 + (qJD(1) * t163 + t161 * t207) * t136;
t103 = -0.2e1 * t162 * t189 * t206 + (t162 * t145 * t186 + (-0.2e1 * t175 * t184 + t180 * t196 + (t124 * t185 + (t123 * t153 - t133 * t180) * t151) * t129) * t144) * t120;
t102 = (t107 * t201 - t115 * t145) * t153 * t178 + ((-t115 * t180 + (-t107 * t147 - t105) * t200) * t145 + (-t115 * t186 - (-t106 * t126 * t152 - t205 * t125 + (t110 * t198 - t125 * t147 - t126 * t179) * t119) * t176 + (t116 * t180 + t153 * t177) * t107 - ((t106 - t179) * t125 + ((-t119 * t152 + 0.1e1) * t147 + (t119 - t152) * t110) * t126) * t145 * t200) * t144) * t111;
t1 = [t141 * t153 * t167 + (t147 * t163 - t180 * t190) * t136, 0, t106, t106, 0, 0; (t115 * t168 + (-t115 * t188 + (qJD(1) * t109 + t105) * t201) * t111) * t152 + (t116 * t168 * t109 + (-((t167 - t188 + (t110 * t141 * t191 + t188) * t136) * t125 + (t166 * t194 - t110 * t144 + (-t139 * t173 + (t110 - 0.2e1 * t187) * t144) * t136) * t126) * t176 + (-t116 * t188 + t144 * t177) * t109 + (-t115 + ((-t148 + t149) * t165 + t204 * t174) * t116) * t144 * qJD(1)) * t111) * t153, 0, t102, t102, 0, 0; 0.2e1 * (-t128 * t131 + t132 * t195) * t206 + ((-t133 * qJD(1) + t150 * t172) * t128 + 0.2e1 * t132 * t175 + (-t131 * t124 - (-t134 * qJD(1) + t151 * t172) * t133 - t132 * t123) * t129) * t120, 0, t103, t103, 0, 0;];
JaD_rot  = t1;
