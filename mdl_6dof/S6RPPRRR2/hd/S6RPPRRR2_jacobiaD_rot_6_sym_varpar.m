% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:35:25
% EndTime: 2019-02-26 20:35:25
% DurationCPUTime: 0.74s
% Computational Cost: add. (3952->98), mult. (2734->203), div. (498->12), fcn. (3199->9), ass. (0->98)
t172 = pkin(11) + qJ(4);
t166 = sin(t172);
t160 = t166 ^ 2;
t168 = cos(t172);
t163 = 0.1e1 / t168 ^ 2;
t221 = t160 * t163;
t174 = qJ(1) + pkin(10);
t167 = sin(t174);
t239 = 0.2e1 * t167;
t238 = t166 * t221;
t169 = cos(t174);
t175 = qJ(5) + qJ(6);
t171 = cos(t175);
t213 = t169 * t171;
t170 = sin(t175);
t216 = t167 * t170;
t149 = t168 * t213 + t216;
t173 = qJD(5) + qJD(6);
t191 = t168 * t173 - qJD(1);
t210 = qJD(4) * t166;
t237 = t191 * t170 + t171 * t210;
t217 = t167 * t166;
t152 = atan2(-t217, -t168);
t151 = cos(t152);
t150 = sin(t152);
t201 = t150 * t217;
t136 = -t151 * t168 - t201;
t133 = 0.1e1 / t136;
t143 = 0.1e1 / t149;
t162 = 0.1e1 / t168;
t134 = 0.1e1 / t136 ^ 2;
t144 = 0.1e1 / t149 ^ 2;
t236 = -0.2e1 * t166;
t161 = t167 ^ 2;
t156 = t161 * t221 + 0.1e1;
t154 = 0.1e1 / t156;
t235 = t154 - 0.1e1;
t211 = qJD(1) * t169;
t198 = t166 * t211;
t208 = qJD(4) * t168;
t209 = qJD(4) * t167;
t127 = (-(-t167 * t208 - t198) * t162 + t209 * t221) * t154;
t223 = t151 * t166;
t122 = (-t127 * t167 + qJD(4)) * t223 + (-t198 + (t127 - t209) * t168) * t150;
t234 = t122 * t133 * t134;
t184 = t168 * t216 + t213;
t197 = t170 * t210;
t128 = t184 * qJD(1) - t149 * t173 + t169 * t197;
t214 = t169 * t170;
t215 = t167 * t171;
t148 = t168 * t214 - t215;
t142 = t148 ^ 2;
t141 = t142 * t144 + 0.1e1;
t225 = t144 * t148;
t190 = -qJD(1) * t168 + t173;
t186 = t190 * t171;
t129 = t167 * t186 - t237 * t169;
t230 = t129 * t143 * t144;
t233 = (-t128 * t225 - t142 * t230) / t141 ^ 2;
t232 = t127 * t150;
t231 = t127 * t166;
t229 = t134 * t166;
t228 = t134 * t169;
t219 = t162 * t166;
t183 = qJD(4) * (t162 * t238 + t219);
t188 = t160 * t167 * t211;
t227 = (t161 * t183 + t163 * t188) / t156 ^ 2;
t195 = 0.1e1 + t221;
t138 = t195 * t167 * t154;
t226 = t138 * t167;
t224 = t150 * t168;
t222 = t160 * t162;
t165 = t169 ^ 2;
t220 = t160 * t165;
t218 = t166 * t169;
t212 = qJD(1) * t167;
t207 = qJD(4) * t169;
t132 = t134 * t220 + 0.1e1;
t206 = 0.2e1 * (-t220 * t234 + (t165 * t166 * t208 - t188) * t134) / t132 ^ 2;
t205 = 0.2e1 * t234;
t204 = -0.2e1 * t233;
t203 = t148 * t230;
t202 = t134 * t218;
t200 = t154 * t222;
t194 = t166 * t206;
t193 = t227 * t236;
t192 = t227 * t239;
t189 = t167 * t200;
t187 = t195 * t169;
t185 = -t143 * t170 + t171 * t225;
t147 = -t168 * t215 + t214;
t139 = 0.1e1 / t141;
t130 = 0.1e1 / t132;
t126 = (t235 * t166 * t150 - t151 * t189) * t169;
t125 = -t167 * t224 + t223 + (-t151 * t217 + t224) * t138;
t123 = -t195 * t192 + (qJD(1) * t187 + t183 * t239) * t154;
t120 = t204 + 0.2e1 * (-t128 * t139 * t144 + (-t139 * t230 - t144 * t233) * t148) * t148;
t1 = [t162 * t169 * t193 + (qJD(4) * t187 - t212 * t219) * t154, 0, 0, t123, 0, 0; (t133 * t194 + (-t133 * t208 + (qJD(1) * t126 + t122) * t229) * t130) * t167 + (t134 * t194 * t126 + (-((t127 * t189 + t235 * t208 + t193) * t150 + (t192 * t222 - t231 + (t231 + (t236 - t238) * t209) * t154) * t151) * t202 + (-t134 * t208 + t166 * t205) * t126 + (-t133 + ((-t161 + t165) * t151 * t200 + t235 * t201) * t134) * t166 * qJD(1)) * t130) * t169, 0, 0 (t125 * t229 - t133 * t168) * t169 * t206 + ((-t133 * t212 + (-qJD(4) * t125 - t122) * t228) * t168 + (-t133 * t207 - (-t123 * t151 * t167 + t150 * t209 + t226 * t232 - t232 + (-qJD(4) * t150 - t151 * t211) * t138) * t202 + (t134 * t212 + t169 * t205) * t125 - ((t123 - t211) * t150 + ((0.1e1 - t226) * qJD(4) + (t138 - t167) * t127) * t151) * t168 * t228) * t166) * t130, 0, 0; 0.2e1 * (t143 * t184 + t147 * t225) * t233 + (0.2e1 * t147 * t203 + (t147 * t128 + t184 * t129 + (-t237 * t167 - t169 * t186) * t148) * t144 + (t190 * t214 + (-t191 * t171 + t197) * t167) * t143) * t139, 0, 0, t185 * t204 * t218 + (t185 * t168 * t207 + (-t185 * t212 + ((-t143 * t173 - 0.2e1 * t203) * t171 + (-t128 * t171 + (-t148 * t173 + t129) * t170) * t144) * t169) * t166) * t139, t120, t120;];
JaD_rot  = t1;
