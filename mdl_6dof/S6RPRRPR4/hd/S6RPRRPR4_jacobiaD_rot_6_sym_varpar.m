% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JaD_rot = S6RPRRPR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:37
% EndTime: 2019-02-26 21:02:38
% DurationCPUTime: 0.76s
% Computational Cost: add. (6073->97), mult. (3810->205), div. (753->12), fcn. (4455->9), ass. (0->96)
t177 = sin(qJ(1));
t238 = 0.2e1 * t177;
t175 = t177 ^ 2;
t172 = pkin(10) + qJ(3) + qJ(4);
t168 = sin(t172);
t164 = t168 ^ 2;
t169 = cos(t172);
t166 = 0.1e1 / t169 ^ 2;
t223 = t164 * t166;
t159 = t175 * t223 + 0.1e1;
t157 = 0.1e1 / t159;
t165 = 0.1e1 / t169;
t178 = cos(qJ(1));
t209 = qJD(1) * t178;
t199 = t168 * t209;
t174 = qJD(3) + qJD(4);
t217 = t174 * t177;
t202 = t166 * t217;
t131 = (-(-t169 * t217 - t199) * t165 + t164 * t202) * t157;
t237 = t131 - t217;
t173 = pkin(11) + qJ(6);
t171 = cos(t173);
t211 = t178 * t171;
t170 = sin(t173);
t214 = t177 * t170;
t153 = t169 * t211 + t214;
t215 = t177 * t168;
t156 = atan2(-t215, -t169);
t155 = cos(t156);
t154 = sin(t156);
t203 = t154 * t215;
t141 = -t155 * t169 - t203;
t138 = 0.1e1 / t141;
t147 = 0.1e1 / t153;
t139 = 0.1e1 / t141 ^ 2;
t148 = 0.1e1 / t153 ^ 2;
t236 = t157 - 0.1e1;
t225 = t155 * t168;
t126 = (-t131 * t177 + t174) * t225 + (t237 * t169 - t199) * t154;
t235 = t126 * t138 * t139;
t188 = t169 * t214 + t211;
t216 = t174 * t178;
t200 = t168 * t216;
t132 = t188 * qJD(1) - t153 * qJD(6) + t170 * t200;
t212 = t178 * t170;
t213 = t177 * t171;
t152 = t169 * t212 - t213;
t146 = t152 ^ 2;
t145 = t146 * t148 + 0.1e1;
t228 = t148 * t152;
t193 = -qJD(1) * t169 + qJD(6);
t194 = qJD(6) * t169 - qJD(1);
t133 = -t194 * t212 + (t177 * t193 - t200) * t171;
t233 = t133 * t147 * t148;
t234 = (-t132 * t228 - t146 * t233) / t145 ^ 2;
t163 = t168 * t164;
t220 = t165 * t168;
t187 = t174 * (t163 * t165 * t166 + t220);
t221 = t164 * t177;
t191 = t209 * t221;
t232 = (t166 * t191 + t175 * t187) / t159 ^ 2;
t231 = t139 * t168;
t230 = t139 * t178;
t229 = t147 * t170;
t227 = t152 * t171;
t226 = t154 * t177;
t224 = t164 * t165;
t176 = t178 ^ 2;
t222 = t164 * t176;
t219 = t168 * t178;
t218 = t169 * t174;
t210 = qJD(1) * t177;
t136 = t139 * t222 + 0.1e1;
t208 = 0.2e1 * (-t222 * t235 + (t168 * t176 * t218 - t191) * t139) / t136 ^ 2;
t207 = 0.2e1 * t235;
t206 = -0.2e1 * t234;
t205 = t139 * t219;
t204 = t152 * t233;
t198 = 0.1e1 + t223;
t197 = t168 * t208;
t196 = -0.2e1 * t168 * t232;
t195 = t232 * t238;
t192 = t155 * t157 * t224;
t190 = t198 * t178;
t189 = t148 * t227 - t229;
t186 = t174 * t215 + t178 * t193;
t151 = -t169 * t213 + t212;
t143 = 0.1e1 / t145;
t142 = t198 * t177 * t157;
t134 = 0.1e1 / t136;
t130 = (t154 * t168 * t236 - t177 * t192) * t178;
t129 = -t169 * t226 + t225 + (t154 * t169 - t155 * t215) * t142;
t127 = -t198 * t195 + (qJD(1) * t190 + t187 * t238) * t157;
t124 = t189 * t206 * t219 + (t189 * t169 * t216 + (-t189 * t210 + ((-qJD(6) * t147 - 0.2e1 * t204) * t171 + (-t132 * t171 + (-qJD(6) * t152 + t133) * t170) * t148) * t178) * t168) * t143;
t123 = (t129 * t231 - t138 * t169) * t178 * t208 + ((-t138 * t210 + (-t129 * t174 - t126) * t230) * t169 + (-t138 * t216 - (-t127 * t155 * t177 - t237 * t154 + (t131 * t226 - t154 * t174 - t155 * t209) * t142) * t205 + (t139 * t210 + t178 * t207) * t129 - ((t127 - t209) * t154 + ((-t142 * t177 + 0.1e1) * t174 + (t142 - t177) * t131) * t155) * t169 * t230) * t168) * t134;
t1 = [t178 * t165 * t196 + (t174 * t190 - t210 * t220) * t157, 0, t127, t127, 0, 0; (t138 * t197 + (-t138 * t218 + (qJD(1) * t130 + t126) * t231) * t134) * t177 + (t139 * t197 * t130 + (-((t196 - t218 + (t131 * t165 * t221 + t218) * t157) * t154 + (t195 * t224 - t131 * t168 + (-t163 * t202 + (t131 - 0.2e1 * t217) * t168) * t157) * t155) * t205 + (-t139 * t218 + t168 * t207) * t130 + (-t138 + ((-t175 + t176) * t192 + t236 * t203) * t139) * t168 * qJD(1)) * t134) * t178, 0, t123, t123, 0, 0; 0.2e1 * (t147 * t188 + t151 * t228) * t234 + (0.2e1 * t151 * t204 - t194 * t147 * t213 + t186 * t229 + (-t152 * t194 * t214 + t151 * t132 + t133 * t188 - t186 * t227) * t148) * t143, 0, t124, t124, 0, t206 + 0.2e1 * (-t132 * t148 * t143 + (-t143 * t233 - t148 * t234) * t152) * t152;];
JaD_rot  = t1;
