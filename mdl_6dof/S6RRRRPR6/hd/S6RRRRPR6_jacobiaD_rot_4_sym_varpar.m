% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR6_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:33:30
% EndTime: 2019-02-26 22:33:31
% DurationCPUTime: 0.76s
% Computational Cost: add. (1570->96), mult. (2734->204), div. (498->12), fcn. (3199->9), ass. (0->95)
t152 = sin(qJ(2));
t145 = t152 ^ 2;
t154 = cos(qJ(2));
t148 = 0.1e1 / t154 ^ 2;
t199 = t145 * t148;
t153 = sin(qJ(1));
t217 = 0.2e1 * t153;
t216 = t152 * t199;
t151 = qJ(3) + qJ(4);
t142 = cos(t151);
t155 = cos(qJ(1));
t191 = t154 * t155;
t141 = sin(t151);
t195 = t153 * t141;
t131 = t142 * t191 + t195;
t193 = t153 * t152;
t135 = atan2(-t193, -t154);
t134 = cos(t135);
t133 = sin(t135);
t180 = t133 * t193;
t121 = -t134 * t154 - t180;
t118 = 0.1e1 / t121;
t125 = 0.1e1 / t131;
t147 = 0.1e1 / t154;
t119 = 0.1e1 / t121 ^ 2;
t126 = 0.1e1 / t131 ^ 2;
t215 = -0.2e1 * t152;
t146 = t153 ^ 2;
t139 = t146 * t199 + 0.1e1;
t137 = 0.1e1 / t139;
t214 = t137 - 0.1e1;
t143 = qJD(3) + qJD(4);
t192 = t153 * t154;
t165 = t141 * t192 + t142 * t155;
t186 = qJD(2) * t155;
t176 = t152 * t186;
t109 = t165 * qJD(1) - t131 * t143 + t141 * t176;
t194 = t153 * t142;
t130 = t141 * t191 - t194;
t124 = t130 ^ 2;
t117 = t124 * t126 + 0.1e1;
t204 = t126 * t130;
t170 = -qJD(1) * t154 + t143;
t171 = t143 * t154 - qJD(1);
t201 = t141 * t155;
t110 = -t171 * t201 + (t170 * t153 - t176) * t142;
t211 = t110 * t125 * t126;
t213 = (-t109 * t204 - t124 * t211) / t117 ^ 2;
t189 = qJD(1) * t155;
t177 = t152 * t189;
t187 = qJD(2) * t154;
t188 = qJD(2) * t153;
t111 = (-(-t153 * t187 - t177) * t147 + t188 * t199) * t137;
t202 = t134 * t152;
t105 = (-t111 * t153 + qJD(2)) * t202 + (-t177 + (t111 - t188) * t154) * t133;
t212 = t105 * t118 * t119;
t210 = t111 * t133;
t209 = t111 * t152;
t208 = t119 * t152;
t197 = t147 * t152;
t164 = qJD(2) * (t147 * t216 + t197);
t168 = t145 * t153 * t189;
t207 = (t146 * t164 + t148 * t168) / t139 ^ 2;
t175 = 0.1e1 + t199;
t123 = t175 * t153 * t137;
t206 = t123 * t153;
t205 = t125 * t141;
t203 = t130 * t142;
t200 = t145 * t147;
t150 = t155 ^ 2;
t198 = t145 * t150;
t196 = t152 * t155;
t190 = qJD(1) * t153;
t114 = t119 * t198 + 0.1e1;
t185 = 0.2e1 * (-t198 * t212 + (t150 * t152 * t187 - t168) * t119) / t114 ^ 2;
t184 = -0.2e1 * t213;
t183 = 0.2e1 * t212;
t182 = t130 * t211;
t181 = t119 * t196;
t179 = t137 * t200;
t174 = t152 * t185;
t173 = t207 * t215;
t172 = t207 * t217;
t169 = t153 * t179;
t167 = t175 * t155;
t166 = t126 * t203 - t205;
t163 = t152 * t188 + t170 * t155;
t129 = -t142 * t192 + t201;
t115 = 0.1e1 / t117;
t112 = 0.1e1 / t114;
t108 = (t214 * t152 * t133 - t134 * t169) * t155;
t107 = -t133 * t192 + t202 + (t133 * t154 - t134 * t193) * t123;
t106 = -t175 * t172 + (qJD(1) * t167 + t164 * t217) * t137;
t102 = t184 + 0.2e1 * (-t109 * t115 * t126 + (-t115 * t211 - t126 * t213) * t130) * t130;
t1 = [t147 * t155 * t173 + (qJD(2) * t167 - t190 * t197) * t137, t106, 0, 0, 0, 0; (t118 * t174 + (-t118 * t187 + (qJD(1) * t108 + t105) * t208) * t112) * t153 + (t119 * t174 * t108 + (-((t111 * t169 + t214 * t187 + t173) * t133 + (t172 * t200 - t209 + (t209 + (t215 - t216) * t188) * t137) * t134) * t181 + (-t119 * t187 + t152 * t183) * t108 + (-t118 + ((-t146 + t150) * t134 * t179 + t214 * t180) * t119) * t152 * qJD(1)) * t112) * t155 (t107 * t208 - t118 * t154) * t155 * t185 + ((-t118 * t190 + (-qJD(2) * t107 - t105) * t155 * t119) * t154 + (-t118 * t186 - (-t106 * t134 * t153 + t133 * t188 + t206 * t210 - t210 + (-qJD(2) * t133 - t134 * t189) * t123) * t181 + (t119 * t190 + t155 * t183) * t107 - ((t106 - t189) * t133 + ((0.1e1 - t206) * qJD(2) + (t123 - t153) * t111) * t134) * t119 * t191) * t152) * t112, 0, 0, 0, 0; 0.2e1 * (t125 * t165 + t129 * t204) * t213 + (0.2e1 * t129 * t182 - t171 * t125 * t194 + t163 * t205 + (-t171 * t130 * t195 + t129 * t109 + t110 * t165 - t163 * t203) * t126) * t115, t166 * t184 * t196 + (t166 * t154 * t186 + (-t166 * t190 + ((-t125 * t143 - 0.2e1 * t182) * t142 + (-t109 * t142 + (-t130 * t143 + t110) * t141) * t126) * t155) * t152) * t115, t102, t102, 0, 0;];
JaD_rot  = t1;
