% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR4_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:52
% EndTime: 2019-02-26 20:05:53
% DurationCPUTime: 0.71s
% Computational Cost: add. (2420->88), mult. (7430->188), div. (524->12), fcn. (9601->11), ass. (0->89)
t164 = sin(pkin(11));
t166 = cos(pkin(11));
t169 = sin(qJ(2));
t167 = cos(pkin(6));
t171 = cos(qJ(2));
t193 = t167 * t171;
t154 = -t164 * t169 + t166 * t193;
t146 = t154 * qJD(2);
t194 = t167 * t169;
t155 = t164 * t171 + t166 * t194;
t168 = sin(qJ(3));
t165 = sin(pkin(6));
t197 = t165 * t168;
t186 = t166 * t197;
t170 = cos(qJ(3));
t191 = qJD(3) * t170;
t121 = -qJD(3) * t186 + t146 * t168 + t155 * t191;
t196 = t165 * t170;
t139 = t155 * t168 + t166 * t196;
t134 = t139 ^ 2;
t158 = -t167 * t170 + t169 * t197;
t152 = 0.1e1 / t158 ^ 2;
t132 = t134 * t152 + 0.1e1;
t129 = 0.1e1 / t132;
t159 = t167 * t168 + t169 * t196;
t195 = t165 * t171;
t185 = qJD(2) * t195;
t144 = t159 * qJD(3) + t168 * t185;
t151 = 0.1e1 / t158;
t202 = t139 * t152;
t109 = (-t121 * t151 + t144 * t202) * t129;
t133 = atan2(-t139, t158);
t125 = sin(t133);
t126 = cos(t133);
t183 = -t125 * t158 - t126 * t139;
t106 = t183 * t109 - t125 * t121 + t126 * t144;
t120 = -t125 * t139 + t126 * t158;
t117 = 0.1e1 / t120;
t118 = 0.1e1 / t120 ^ 2;
t214 = t106 * t117 * t118;
t179 = t164 * t194 - t166 * t171;
t181 = t164 * t196 + t168 * t179;
t213 = -0.2e1 * t181;
t212 = 0.2e1 * t170;
t211 = t213 * t214;
t178 = -t151 * t154 + t195 * t202;
t210 = t168 * t178;
t201 = t144 * t151 * t152;
t209 = -0.2e1 * (t121 * t202 - t134 * t201) / t132 ^ 2;
t143 = t164 * t197 - t170 * t179;
t136 = 0.1e1 / t143;
t137 = 0.1e1 / t143 ^ 2;
t180 = t164 * t193 + t166 * t169;
t150 = t180 ^ 2;
t199 = t150 * t137;
t131 = 0.1e1 + t199;
t148 = t180 * qJD(2);
t124 = t181 * qJD(3) - t148 * t170;
t205 = t124 * t136 * t137;
t187 = t150 * t205;
t149 = t179 * qJD(2);
t200 = t149 * t180;
t208 = (-t137 * t200 - t187) / t131 ^ 2;
t207 = t118 * t181;
t123 = t143 * qJD(3) - t148 * t168;
t206 = t123 * t118;
t204 = t125 * t181;
t203 = t126 * t181;
t198 = t180 * t168;
t192 = qJD(2) * t169;
t135 = t181 ^ 2;
t115 = t135 * t118 + 0.1e1;
t190 = 0.2e1 * (-t135 * t214 - t181 * t206) / t115 ^ 2;
t188 = t180 * t213;
t184 = -0.2e1 * t139 * t201;
t141 = t155 * t170 - t186;
t182 = -t141 * t151 + t159 * t202;
t147 = t155 * qJD(2);
t145 = -t158 * qJD(3) + t170 * t185;
t127 = 0.1e1 / t131;
t122 = -t139 * qJD(3) + t146 * t170;
t112 = 0.1e1 / t115;
t111 = t129 * t210;
t110 = t182 * t129;
t108 = (-t125 * t154 + t126 * t195) * t168 + t183 * t111;
t107 = t183 * t110 - t125 * t141 + t126 * t159;
t105 = t182 * t209 + (t159 * t184 - t122 * t151 + (t121 * t159 + t139 * t145 + t141 * t144) * t152) * t129;
t103 = t209 * t210 + (t178 * t191 + (t184 * t195 + t147 * t151 + (t144 * t154 + (t121 * t171 - t139 * t192) * t165) * t152) * t168) * t129;
t1 = [0, t103, t105, 0, 0, 0; 0 (-t108 * t207 + t117 * t198) * t190 + ((t149 * t168 - t180 * t191) * t117 + (-t206 + t211) * t108 + (t198 * t106 + (-t103 * t139 - t111 * t121 + (-t168 * t192 + t171 * t191) * t165 + (-t111 * t158 - t154 * t168) * t109) * t203 + (-t154 * t191 - t103 * t158 - t111 * t144 + t147 * t168 + (t111 * t139 - t168 * t195) * t109) * t204) * t118) * t112 (-t107 * t207 - t117 * t143) * t190 + (t107 * t211 + t124 * t117 + (-t143 * t106 - t107 * t123 + (-t105 * t139 - t110 * t121 + t145 + (-t110 * t158 - t141) * t109) * t203 + (-t105 * t158 - t110 * t144 - t122 + (t110 * t139 - t159) * t109) * t204) * t118) * t112, 0, 0, 0; 0, 0.2e1 * (-t136 * t179 + t170 * t199) * t208 + (t187 * t212 + t136 * t148 + (qJD(3) * t150 * t168 - t124 * t179 + t200 * t212) * t137) * t127, t137 * t188 * t208 + (t188 * t205 + (-t123 * t180 - t149 * t181) * t137) * t127, 0, 0, 0;];
JaD_rot  = t1;
