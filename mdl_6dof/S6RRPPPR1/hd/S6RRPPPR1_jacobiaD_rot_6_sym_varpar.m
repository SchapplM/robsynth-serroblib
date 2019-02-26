% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:49
% EndTime: 2019-02-26 21:21:50
% DurationCPUTime: 0.83s
% Computational Cost: add. (2357->105), mult. (3345->232), div. (468->12), fcn. (4023->11), ass. (0->107)
t168 = qJ(2) + pkin(9);
t166 = sin(t168);
t162 = t166 ^ 2;
t167 = cos(t168);
t164 = 0.1e1 / t167 ^ 2;
t223 = t162 * t164;
t174 = sin(qJ(1));
t169 = t174 ^ 2;
t160 = t169 * t223 + 0.1e1;
t163 = 0.1e1 / t167;
t220 = t163 * t166;
t242 = t166 * t223;
t184 = qJD(2) * (t163 * t242 + t220);
t176 = cos(qJ(1));
t212 = qJD(1) * t176;
t221 = t162 * t174;
t190 = t212 * t221;
t228 = (t164 * t190 + t169 * t184) / t160 ^ 2;
t243 = -0.2e1 * t228;
t194 = 0.1e1 + t223;
t241 = t174 * t194;
t157 = 0.1e1 / t160;
t201 = t166 * t212;
t210 = qJD(2) * t174;
t126 = ((t167 * t210 + t201) * t163 + t210 * t223) * t157;
t240 = t126 - t210;
t173 = sin(qJ(6));
t175 = cos(qJ(6));
t208 = qJD(6) * t176;
t213 = qJD(1) * t174;
t239 = t173 * t213 - t175 * t208;
t238 = t173 * t208 + t175 * t213;
t216 = t174 * t166;
t159 = atan2(t216, t167);
t156 = cos(t159);
t155 = sin(t159);
t203 = t155 * t216;
t133 = t156 * t167 + t203;
t130 = 0.1e1 / t133;
t172 = cos(pkin(10));
t214 = t174 * t172;
t171 = sin(pkin(10));
t218 = t171 * t176;
t153 = t167 * t218 - t214;
t215 = t174 * t171;
t217 = t172 * t176;
t154 = t167 * t217 + t215;
t143 = t153 * t173 + t154 * t175;
t137 = 0.1e1 / t143;
t131 = 0.1e1 / t133 ^ 2;
t138 = 0.1e1 / t143 ^ 2;
t237 = 0.2e1 * t166;
t236 = t157 - 0.1e1;
t170 = t176 ^ 2;
t222 = t162 * t170;
t129 = t131 * t222 + 0.1e1;
t211 = qJD(2) * t167;
t224 = t156 * t166;
t117 = (t126 * t174 - qJD(2)) * t224 + (-t240 * t167 + t201) * t155;
t233 = t117 * t130 * t131;
t235 = (-t222 * t233 + (t166 * t170 * t211 - t190) * t131) / t129 ^ 2;
t151 = -t167 * t215 - t217;
t209 = qJD(2) * t176;
t197 = t166 * t209;
t144 = qJD(1) * t151 - t171 * t197;
t152 = -t167 * t214 + t218;
t145 = qJD(1) * t152 - t172 * t197;
t120 = qJD(6) * t143 - t144 * t175 + t145 * t173;
t187 = t153 * t175 - t154 * t173;
t136 = t187 ^ 2;
t125 = t136 * t138 + 0.1e1;
t227 = t138 * t187;
t121 = qJD(6) * t187 + t144 * t173 + t145 * t175;
t139 = t137 * t138;
t232 = t121 * t139;
t234 = (-t120 * t227 - t136 * t232) / t125 ^ 2;
t231 = t126 * t166;
t230 = t131 * t166;
t229 = t131 * t176;
t185 = -t171 * t173 - t172 * t175;
t219 = t166 * t176;
t149 = t185 * t219;
t226 = t138 * t149;
t225 = t155 * t174;
t207 = 0.2e1 * t234;
t206 = -0.2e1 * t233;
t205 = -0.2e1 * t139 * t187;
t204 = t131 * t219;
t202 = t157 * t162 * t163;
t198 = t166 * t210;
t193 = -0.2e1 * t166 * t235;
t192 = t163 * t243;
t191 = t174 * t202;
t189 = t194 * t176;
t188 = t151 * t175 - t152 * t173;
t141 = t151 * t173 + t152 * t175;
t186 = t171 * t175 - t172 * t173;
t148 = t186 * t219;
t147 = -qJD(1) * t154 + t172 * t198;
t146 = -qJD(1) * t153 + t171 * t198;
t135 = t157 * t241;
t127 = 0.1e1 / t129;
t123 = 0.1e1 / t125;
t122 = (-t155 * t166 * t236 + t156 * t191) * t176;
t119 = t167 * t225 - t224 + (-t155 * t167 + t156 * t216) * t135;
t118 = t241 * t243 + (qJD(1) * t189 + 0.2e1 * t174 * t184) * t157;
t1 = [t192 * t219 + (qJD(2) * t189 - t213 * t220) * t157, t118, 0, 0, 0, 0; (t130 * t193 + (t130 * t211 + (-qJD(1) * t122 - t117) * t230) * t127) * t174 + (t131 * t193 * t122 + (((-t126 * t191 - t211 * t236 + t228 * t237) * t155 + (t192 * t221 + t231 + (-t231 + (t237 + t242) * t210) * t157) * t156) * t204 + (t131 * t211 + t166 * t206) * t122 + (t130 + ((-t169 + t170) * t156 * t202 + t236 * t203) * t131) * t166 * qJD(1)) * t127) * t176, 0.2e1 * (-t119 * t230 + t130 * t167) * t176 * t235 + ((t130 * t213 + (qJD(2) * t119 + t117) * t229) * t167 + (t130 * t209 + (t118 * t156 * t174 + t240 * t155 + (qJD(2) * t155 - t126 * t225 + t156 * t212) * t135) * t204 + (-t131 * t213 + t176 * t206) * t119 + ((-t118 + t212) * t155 + ((t135 * t174 - 0.1e1) * qJD(2) + (-t135 + t174) * t126) * t156) * t167 * t229) * t166) * t127, 0, 0, 0, 0; (t137 * t188 - t141 * t227) * t207 + ((qJD(6) * t141 - t146 * t175 + t147 * t173) * t137 + t141 * t121 * t205 + (t188 * t121 + (qJD(6) * t188 + t146 * t173 + t147 * t175) * t187 - t141 * t120) * t138) * t123 (-t137 * t148 - t187 * t226) * t207 + (-t120 * t226 + (-t138 * t148 + t149 * t205) * t121 + (t137 * t186 + t185 * t227) * t167 * t209 + ((t239 * t137 + t238 * t227) * t172 + (-t238 * t137 + t239 * t227) * t171) * t166) * t123, 0, 0, 0, -0.2e1 * t234 - 0.2e1 * (t120 * t123 * t138 - (-t123 * t232 - t138 * t234) * t187) * t187;];
JaD_rot  = t1;
