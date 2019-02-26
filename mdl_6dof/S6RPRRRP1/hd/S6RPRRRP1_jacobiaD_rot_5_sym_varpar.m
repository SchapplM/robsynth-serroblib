% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:08:01
% EndTime: 2019-02-26 21:08:02
% DurationCPUTime: 0.81s
% Computational Cost: add. (5102->98), mult. (3810->203), div. (753->12), fcn. (4455->9), ass. (0->97)
t173 = qJ(3) + qJ(4);
t169 = sin(t173);
t163 = t169 ^ 2;
t170 = cos(t173);
t165 = 0.1e1 / t170 ^ 2;
t220 = t163 * t165;
t172 = qJ(1) + pkin(10);
t167 = sin(t172);
t239 = 0.2e1 * t167;
t238 = t169 * t220;
t168 = cos(t172);
t174 = sin(qJ(5));
t175 = cos(qJ(5));
t211 = t170 * t175;
t150 = t167 * t174 + t168 * t211;
t192 = qJD(5) * t170 - qJD(1);
t171 = qJD(3) + qJD(4);
t214 = t169 * t171;
t237 = t192 * t174 + t175 * t214;
t218 = t167 * t169;
t153 = atan2(-t218, -t170);
t152 = cos(t153);
t151 = sin(t153);
t201 = t151 * t218;
t138 = -t152 * t170 - t201;
t135 = 0.1e1 / t138;
t144 = 0.1e1 / t150;
t164 = 0.1e1 / t170;
t136 = 0.1e1 / t138 ^ 2;
t145 = 0.1e1 / t150 ^ 2;
t236 = -0.2e1 * t169;
t160 = t167 ^ 2;
t156 = t160 * t220 + 0.1e1;
t154 = 0.1e1 / t156;
t235 = t154 - 0.1e1;
t207 = qJD(1) * t169;
t197 = t168 * t207;
t213 = t170 * t171;
t217 = t167 * t171;
t128 = (-(-t167 * t213 - t197) * t164 + t217 * t220) * t154;
t223 = t152 * t169;
t123 = (-t128 * t167 + t171) * t223 + (-t197 + (t128 - t217) * t170) * t151;
t234 = t123 * t135 * t136;
t233 = t128 * t151;
t232 = t128 * t169;
t191 = -qJD(1) * t170 + qJD(5);
t187 = t191 * t175;
t133 = t167 * t187 - t237 * t168;
t231 = t133 * t144 * t145;
t184 = (t169 + t238) * t164 * t171;
t208 = qJD(1) * t168;
t189 = t163 * t167 * t208;
t230 = (t160 * t184 + t165 * t189) / t156 ^ 2;
t229 = t136 * t168;
t228 = t136 * t169;
t196 = 0.1e1 + t220;
t139 = t196 * t167 * t154;
t227 = t139 * t167;
t212 = t170 * t174;
t185 = t167 * t212 + t168 * t175;
t199 = t174 * t214;
t132 = t185 * qJD(1) - t150 * qJD(5) + t168 * t199;
t149 = -t167 * t175 + t168 * t212;
t143 = t149 ^ 2;
t142 = t143 * t145 + 0.1e1;
t225 = t145 * t149;
t226 = 0.1e1 / t142 ^ 2 * (-t132 * t225 - t143 * t231);
t224 = t151 * t170;
t161 = t168 ^ 2;
t222 = t161 * t163;
t221 = t163 * t164;
t219 = t164 * t167;
t215 = t168 * t174;
t210 = t171 * t135;
t209 = qJD(1) * t167;
t131 = t136 * t222 + 0.1e1;
t206 = 0.2e1 * (-t222 * t234 + (t161 * t169 * t213 - t189) * t136) / t131 ^ 2;
t205 = 0.2e1 * t234;
t204 = 0.2e1 * t226;
t203 = t149 * t231;
t202 = t168 * t228;
t195 = t169 * t206;
t194 = t230 * t239;
t193 = t230 * t236;
t190 = t152 * t154 * t221;
t188 = t196 * t168;
t186 = -t144 * t174 + t175 * t225;
t183 = t186 * t169;
t148 = -t167 * t211 + t215;
t140 = 0.1e1 / t142;
t129 = 0.1e1 / t131;
t127 = (t235 * t169 * t151 - t167 * t190) * t168;
t126 = -t167 * t224 + t223 + (-t152 * t218 + t224) * t139;
t124 = -t196 * t194 + (qJD(1) * t188 + t184 * t239) * t154;
t121 = -t168 * t183 * t204 + (-t183 * t209 + (t186 * t213 + ((-qJD(5) * t144 - 0.2e1 * t203) * t175 + (-t132 * t175 + (-qJD(5) * t149 + t133) * t174) * t145) * t169) * t168) * t140;
t120 = (t126 * t228 - t135 * t170) * t168 * t206 + ((-t135 * t209 + (-t126 * t171 - t123) * t229) * t170 + (-t168 * t210 - (-t124 * t152 * t167 + t151 * t217 + t227 * t233 - t233 + (-t151 * t171 - t152 * t208) * t139) * t202 + (t136 * t209 + t168 * t205) * t126 - ((t124 - t208) * t151 + ((0.1e1 - t227) * t171 + (t139 - t167) * t128) * t152) * t170 * t229) * t169) * t129;
t1 = [t168 * t164 * t193 + (t171 * t188 - t207 * t219) * t154, 0, t124, t124, 0, 0; (t135 * t195 + (-t170 * t210 + (qJD(1) * t127 + t123) * t228) * t129) * t167 + (t136 * t195 * t127 + (-((t193 - t213 + (t128 * t163 * t219 + t213) * t154) * t151 + (t194 * t221 - t232 + (t232 + (t236 - t238) * t217) * t154) * t152) * t202 + (-t136 * t213 + t169 * t205) * t127 + (-t135 + ((-t160 + t161) * t190 + t235 * t201) * t136) * t207) * t129) * t168, 0, t120, t120, 0, 0; (t144 * t185 + t148 * t225) * t204 + (0.2e1 * t148 * t203 + (t148 * t132 + t185 * t133 + (-t237 * t167 - t168 * t187) * t149) * t145 + (t191 * t215 + (-t192 * t175 + t199) * t167) * t144) * t140, 0, t121, t121, -0.2e1 * t226 + 0.2e1 * (-t132 * t145 * t140 + (-t140 * t231 - t145 * t226) * t149) * t149, 0;];
JaD_rot  = t1;
