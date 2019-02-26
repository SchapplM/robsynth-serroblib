% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRRRRR5
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRRR5_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_jacobiaD_rot_3_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:11
% EndTime: 2019-02-26 20:21:11
% DurationCPUTime: 0.52s
% Computational Cost: add. (1326->63), mult. (4214->161), div. (275->12), fcn. (5416->13), ass. (0->79)
t164 = sin(pkin(7));
t162 = t164 ^ 2;
t206 = 0.2e1 * t162;
t166 = cos(pkin(13));
t163 = sin(pkin(13));
t170 = sin(qJ(2));
t168 = cos(pkin(6));
t172 = cos(qJ(2));
t190 = t168 * t172;
t182 = -t163 * t170 + t166 * t190;
t165 = sin(pkin(6));
t167 = cos(pkin(7));
t195 = t165 * t167;
t146 = t182 * t164 + t166 * t195;
t196 = t164 * t172;
t156 = -t165 * t196 + t168 * t167;
t141 = atan2(t146, t156);
t136 = sin(t141);
t137 = cos(t141);
t123 = t136 * t146 + t137 * t156;
t120 = 0.1e1 / t123;
t169 = sin(qJ(3));
t171 = cos(qJ(3));
t191 = t168 * t170;
t180 = t163 * t191 - t166 * t172;
t181 = t163 * t190 + t166 * t170;
t187 = t163 * t164 * t165;
t183 = -t167 * t181 + t187;
t135 = t183 * t169 - t171 * t180;
t131 = 0.1e1 / t135;
t153 = 0.1e1 / t156;
t121 = 0.1e1 / t123 ^ 2;
t132 = 0.1e1 / t135 ^ 2;
t154 = 0.1e1 / t156 ^ 2;
t147 = t163 * t195 + t164 * t181;
t145 = t147 ^ 2;
t119 = t145 * t121 + 0.1e1;
t152 = t180 * qJD(2);
t157 = -t163 * t172 - t166 * t191;
t150 = t157 * qJD(2);
t194 = t165 * t170;
t198 = t146 * t154;
t185 = t194 * t198;
t144 = t146 ^ 2;
t140 = t144 * t154 + 0.1e1;
t138 = 0.1e1 / t140;
t199 = t138 * t164;
t115 = (-qJD(2) * t185 + t150 * t153) * t199;
t184 = -t136 * t156 + t137 * t146;
t189 = qJD(2) * t165;
t186 = t170 * t189;
t112 = (t136 * t150 + t137 * t186) * t164 + t184 * t115;
t204 = t112 * t120 * t121;
t205 = (-t147 * t121 * t152 * t164 - t145 * t204) / t119 ^ 2;
t192 = t167 * t171;
t197 = t180 * t169;
t134 = -t171 * t187 + t181 * t192 - t197;
t130 = t134 ^ 2;
t127 = t130 * t132 + 0.1e1;
t151 = t181 * qJD(2);
t193 = t167 * t169;
t129 = t152 * t193 - t151 * t171 + (t183 * t171 + t197) * qJD(3);
t201 = t129 * t131 * t132;
t128 = t135 * qJD(3) - t151 * t169 - t152 * t192;
t202 = t128 * t132;
t203 = (-t130 * t201 + t134 * t202) / t127 ^ 2;
t143 = -t171 * t181 + t180 * t193;
t200 = t134 * t143;
t188 = t154 * t162 * t170;
t142 = -t169 * t181 - t180 * t192;
t179 = -t153 * t157 + t185;
t155 = t153 * t154;
t149 = t182 * qJD(2);
t125 = 0.1e1 / t127;
t117 = 0.1e1 / t119;
t116 = t179 * t199;
t113 = (t136 * t157 + t137 * t194) * t164 - t184 * t116;
t111 = t179 / t140 ^ 2 * (-t144 * t155 * t186 + t150 * t198) * t206 + (-t149 * t153 * t164 + (-t150 * t188 + (-t157 * t188 + (t155 * t165 * t170 ^ 2 * t206 - t154 * t196) * t146) * qJD(2)) * t165) * t138;
t1 = [0, t111, 0, 0, 0, 0; 0 (-(t123 * t116 * t115 + t184 * t111) * t121 * t117 + 0.2e1 * (t117 * t204 + t121 * t205) * t113) * t147 + (0.2e1 * t180 * t120 * t205 + (-t151 * t120 + (t180 * t112 + t113 * t152 + (-(t115 * t157 - t116 * t150 + t172 * t189) * t137 - (-t149 + (qJD(2) * t116 - t115) * t194) * t136) * t147) * t121) * t117) * t164, 0, 0, 0, 0; 0, 0.2e1 * (-t131 * t142 + t132 * t200) * t203 + ((t143 * qJD(3) - t151 * t192 + t152 * t169) * t131 + 0.2e1 * t200 * t201 + (-t142 * t129 - (-t142 * qJD(3) + t151 * t193 + t152 * t171) * t134 - t143 * t128) * t132) * t125, -0.2e1 * t203 + 0.2e1 * (t125 * t202 + (-t125 * t201 - t132 * t203) * t134) * t134, 0, 0, 0;];
JaD_rot  = t1;
