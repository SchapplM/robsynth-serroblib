% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:05
% EndTime: 2019-02-26 22:11:06
% DurationCPUTime: 0.73s
% Computational Cost: add. (3360->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->95)
t173 = sin(qJ(1));
t171 = qJ(2) + qJ(3);
t166 = sin(t171);
t162 = 0.1e1 / t166 ^ 2;
t167 = cos(t171);
t165 = t167 ^ 2;
t219 = t162 * t165;
t194 = 0.1e1 + t219;
t234 = t173 * t194;
t169 = t173 ^ 2;
t159 = t169 * t219 + 0.1e1;
t157 = 0.1e1 / t159;
t161 = 0.1e1 / t166;
t175 = cos(qJ(1));
t206 = qJD(1) * t175;
t195 = t167 * t206;
t168 = qJD(2) + qJD(3);
t214 = t168 * t173;
t197 = t162 * t214;
t131 = ((t166 * t214 - t195) * t161 + t165 * t197) * t157;
t233 = -t131 + t214;
t190 = qJD(1) * t166 + qJD(5);
t213 = t168 * t175;
t232 = -t167 * t213 + t190 * t173;
t212 = t173 * t167;
t156 = atan2(-t212, t166);
t155 = cos(t156);
t154 = sin(t156);
t199 = t154 * t212;
t141 = t155 * t166 - t199;
t138 = 0.1e1 / t141;
t172 = sin(qJ(5));
t209 = t175 * t172;
t174 = cos(qJ(5));
t210 = t173 * t174;
t151 = t166 * t209 + t210;
t147 = 0.1e1 / t151;
t139 = 0.1e1 / t141 ^ 2;
t148 = 0.1e1 / t151 ^ 2;
t231 = t157 - 0.1e1;
t221 = t155 * t167;
t126 = (-t131 * t173 + t168) * t221 + (t233 * t166 - t195) * t154;
t230 = t126 * t138 * t139;
t191 = qJD(5) * t166 + qJD(1);
t186 = t191 * t175;
t135 = t172 * t186 + t232 * t174;
t208 = t175 * t174;
t211 = t173 * t172;
t150 = -t166 * t208 + t211;
t146 = t150 ^ 2;
t145 = t146 * t148 + 0.1e1;
t224 = t148 * t150;
t136 = -t232 * t172 + t174 * t186;
t228 = t136 * t147 * t148;
t229 = (t135 * t224 - t146 * t228) / t145 ^ 2;
t164 = t167 * t165;
t220 = t161 * t167;
t184 = t168 * (-t161 * t162 * t164 - t220);
t217 = t165 * t173;
t188 = t206 * t217;
t227 = (t162 * t188 + t169 * t184) / t159 ^ 2;
t226 = t139 * t167;
t225 = t139 * t175;
t223 = t150 * t172;
t222 = t154 * t173;
t170 = t175 ^ 2;
t218 = t165 * t170;
t216 = t166 * t168;
t215 = t167 * t168;
t207 = qJD(1) * t173;
t134 = t139 * t218 + 0.1e1;
t205 = 0.2e1 * (-t218 * t230 + (-t166 * t170 * t215 - t188) * t139) / t134 ^ 2;
t204 = 0.2e1 * t230;
t203 = 0.2e1 * t229;
t202 = -0.2e1 * t227;
t201 = t167 * t227;
t200 = t167 * t225;
t198 = t161 * t217;
t193 = t167 * t205;
t192 = 0.2e1 * t150 * t228;
t189 = t155 * t157 * t161 * t165;
t187 = t194 * t175;
t185 = t147 * t174 + t148 * t223;
t183 = t185 * t175;
t153 = -t166 * t211 + t208;
t152 = t166 * t210 + t209;
t143 = 0.1e1 / t145;
t142 = t157 * t234;
t132 = 0.1e1 / t134;
t130 = (t231 * t167 * t154 + t173 * t189) * t175;
t128 = t166 * t222 + t221 + (-t154 * t166 - t155 * t212) * t142;
t127 = t202 * t234 + (qJD(1) * t187 + 0.2e1 * t173 * t184) * t157;
t124 = t167 * t183 * t203 + (t183 * t216 + (t185 * t207 + ((qJD(5) * t147 + t192) * t172 + (-t135 * t172 + (-qJD(5) * t150 + t136) * t174) * t148) * t175) * t167) * t143;
t123 = (t128 * t226 + t138 * t166) * t175 * t205 + ((t138 * t207 + (t128 * t168 + t126) * t225) * t166 + (-t138 * t213 - (-t127 * t155 * t173 + t233 * t154 + (t131 * t222 - t154 * t168 - t155 * t206) * t142) * t200 + (t139 * t207 + t175 * t204) * t128 - ((-t127 + t206) * t154 + ((t142 * t173 - 0.1e1) * t168 + (-t142 + t173) * t131) * t155) * t166 * t225) * t167) * t132;
t1 = [0.2e1 * t175 * t161 * t201 + (t168 * t187 + t207 * t220) * t157, t127, t127, 0, 0, 0; (t138 * t193 + (t138 * t216 + (qJD(1) * t130 + t126) * t226) * t132) * t173 + (t139 * t193 * t130 + (-((-0.2e1 * t201 + t216 + (-t131 * t198 - t216) * t157) * t154 + (t198 * t202 - t131 * t167 + (-t164 * t197 + (t131 - 0.2e1 * t214) * t167) * t157) * t155) * t200 + (t139 * t216 + t167 * t204) * t130 + (-t138 + ((t169 - t170) * t189 + t231 * t199) * t139) * t167 * qJD(1)) * t132) * t175, t123, t123, 0, 0, 0; (-t147 * t152 + t153 * t224) * t203 + (t153 * t192 + (-t153 * t135 - t152 * t136 + t191 * t150 * t210 - (-t168 * t212 - t190 * t175) * t223) * t148 + (t190 * t208 + (-t191 * t172 + t174 * t215) * t173) * t147) * t143, t124, t124, 0, -0.2e1 * t229 + 0.2e1 * (t135 * t148 * t143 + (-t143 * t228 - t148 * t229) * t150) * t150, 0;];
JaD_rot  = t1;
