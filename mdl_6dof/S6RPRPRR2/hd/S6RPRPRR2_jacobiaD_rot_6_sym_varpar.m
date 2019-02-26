% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:49:39
% EndTime: 2019-02-26 20:49:39
% DurationCPUTime: 0.76s
% Computational Cost: add. (3952->98), mult. (2734->203), div. (498->12), fcn. (3199->9), ass. (0->98)
t176 = qJ(3) + pkin(11);
t169 = sin(t176);
t163 = t169 ^ 2;
t171 = cos(t176);
t166 = 0.1e1 / t171 ^ 2;
t224 = t163 * t166;
t177 = qJ(1) + pkin(10);
t170 = sin(t177);
t242 = 0.2e1 * t170;
t241 = t169 * t224;
t172 = cos(t177);
t178 = qJ(5) + qJ(6);
t174 = cos(t178);
t216 = t172 * t174;
t173 = sin(t178);
t219 = t170 * t173;
t152 = t171 * t216 + t219;
t175 = qJD(5) + qJD(6);
t194 = t171 * t175 - qJD(1);
t213 = qJD(3) * t169;
t240 = t194 * t173 + t174 * t213;
t220 = t170 * t169;
t155 = atan2(-t220, -t171);
t154 = cos(t155);
t153 = sin(t155);
t204 = t153 * t220;
t139 = -t154 * t171 - t204;
t136 = 0.1e1 / t139;
t146 = 0.1e1 / t152;
t165 = 0.1e1 / t171;
t137 = 0.1e1 / t139 ^ 2;
t147 = 0.1e1 / t152 ^ 2;
t239 = -0.2e1 * t169;
t164 = t170 ^ 2;
t159 = t164 * t224 + 0.1e1;
t157 = 0.1e1 / t159;
t238 = t157 - 0.1e1;
t214 = qJD(1) * t172;
t201 = t169 * t214;
t211 = qJD(3) * t171;
t212 = qJD(3) * t170;
t130 = (-(-t170 * t211 - t201) * t165 + t212 * t224) * t157;
t226 = t154 * t169;
t125 = (-t130 * t170 + qJD(3)) * t226 + (-t201 + (t130 - t212) * t171) * t153;
t237 = t125 * t136 * t137;
t187 = t171 * t219 + t216;
t200 = t173 * t213;
t131 = t187 * qJD(1) - t152 * t175 + t172 * t200;
t217 = t172 * t173;
t218 = t170 * t174;
t151 = t171 * t217 - t218;
t145 = t151 ^ 2;
t144 = t145 * t147 + 0.1e1;
t228 = t147 * t151;
t193 = -qJD(1) * t171 + t175;
t189 = t193 * t174;
t132 = t170 * t189 - t240 * t172;
t233 = t132 * t146 * t147;
t236 = (-t131 * t228 - t145 * t233) / t144 ^ 2;
t235 = t130 * t153;
t234 = t130 * t169;
t232 = t137 * t169;
t231 = t137 * t172;
t222 = t165 * t169;
t186 = qJD(3) * (t165 * t241 + t222);
t191 = t163 * t170 * t214;
t230 = (t164 * t186 + t166 * t191) / t159 ^ 2;
t198 = 0.1e1 + t224;
t141 = t198 * t170 * t157;
t229 = t141 * t170;
t227 = t153 * t171;
t225 = t163 * t165;
t168 = t172 ^ 2;
t223 = t163 * t168;
t221 = t169 * t172;
t215 = qJD(1) * t170;
t210 = qJD(3) * t172;
t135 = t137 * t223 + 0.1e1;
t209 = 0.2e1 * (-t223 * t237 + (t168 * t169 * t211 - t191) * t137) / t135 ^ 2;
t208 = 0.2e1 * t237;
t207 = -0.2e1 * t236;
t206 = t151 * t233;
t205 = t137 * t221;
t203 = t157 * t225;
t197 = t169 * t209;
t196 = t230 * t239;
t195 = t230 * t242;
t192 = t170 * t203;
t190 = t198 * t172;
t188 = -t146 * t173 + t174 * t228;
t150 = -t171 * t218 + t217;
t142 = 0.1e1 / t144;
t133 = 0.1e1 / t135;
t129 = (t238 * t169 * t153 - t154 * t192) * t172;
t128 = -t170 * t227 + t226 + (-t154 * t220 + t227) * t141;
t126 = -t198 * t195 + (qJD(1) * t190 + t186 * t242) * t157;
t123 = t207 + 0.2e1 * (-t131 * t142 * t147 + (-t142 * t233 - t147 * t236) * t151) * t151;
t1 = [t165 * t172 * t196 + (qJD(3) * t190 - t215 * t222) * t157, 0, t126, 0, 0, 0; (t136 * t197 + (-t136 * t211 + (qJD(1) * t129 + t125) * t232) * t133) * t170 + (t137 * t197 * t129 + (-((t130 * t192 + t238 * t211 + t196) * t153 + (t195 * t225 - t234 + (t234 + (t239 - t241) * t212) * t157) * t154) * t205 + (-t137 * t211 + t169 * t208) * t129 + (-t136 + ((-t164 + t168) * t154 * t203 + t238 * t204) * t137) * t169 * qJD(1)) * t133) * t172, 0 (t128 * t232 - t136 * t171) * t172 * t209 + ((-t136 * t215 + (-qJD(3) * t128 - t125) * t231) * t171 + (-t136 * t210 - (-t126 * t154 * t170 + t153 * t212 + t229 * t235 - t235 + (-qJD(3) * t153 - t154 * t214) * t141) * t205 + (t137 * t215 + t172 * t208) * t128 - ((t126 - t214) * t153 + ((0.1e1 - t229) * qJD(3) + (t141 - t170) * t130) * t154) * t171 * t231) * t169) * t133, 0, 0, 0; 0.2e1 * (t146 * t187 + t150 * t228) * t236 + (0.2e1 * t150 * t206 + (t150 * t131 + t187 * t132 + (-t240 * t170 - t172 * t189) * t151) * t147 + (t193 * t217 + (-t194 * t174 + t200) * t170) * t146) * t142, 0, t188 * t207 * t221 + (t188 * t171 * t210 + (-t188 * t215 + ((-t146 * t175 - 0.2e1 * t206) * t174 + (-t131 * t174 + (-t151 * t175 + t132) * t173) * t147) * t172) * t169) * t142, 0, t123, t123;];
JaD_rot  = t1;
