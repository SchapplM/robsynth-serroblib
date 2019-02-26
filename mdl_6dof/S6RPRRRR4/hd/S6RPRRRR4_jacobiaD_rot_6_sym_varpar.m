% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:26
% EndTime: 2019-02-26 21:16:27
% DurationCPUTime: 0.81s
% Computational Cost: add. (10386->97), mult. (5101->208), div. (1026->12), fcn. (5942->9), ass. (0->95)
t180 = sin(qJ(1));
t242 = 0.2e1 * t180;
t177 = t180 ^ 2;
t174 = pkin(11) + qJ(3) + qJ(4) + qJ(5);
t171 = sin(t174);
t167 = t171 ^ 2;
t172 = cos(t174);
t169 = 0.1e1 / t172 ^ 2;
t227 = t167 * t169;
t164 = t177 * t227 + 0.1e1;
t158 = 0.1e1 / t164;
t168 = 0.1e1 / t172;
t182 = cos(qJ(1));
t213 = qJD(1) * t182;
t203 = t171 * t213;
t176 = qJD(3) + qJD(4) + qJD(5);
t221 = t176 * t180;
t206 = t169 * t221;
t136 = (-(-t172 * t221 - t203) * t168 + t167 * t206) * t158;
t241 = t136 - t221;
t181 = cos(qJ(6));
t215 = t181 * t182;
t179 = sin(qJ(6));
t217 = t180 * t179;
t163 = t172 * t215 + t217;
t218 = t180 * t171;
t153 = atan2(-t218, -t172);
t152 = cos(t153);
t151 = sin(t153);
t207 = t151 * t218;
t144 = -t152 * t172 - t207;
t141 = 0.1e1 / t144;
t155 = 0.1e1 / t163;
t142 = 0.1e1 / t144 ^ 2;
t156 = 0.1e1 / t163 ^ 2;
t240 = t158 - 0.1e1;
t232 = t152 * t171;
t131 = (-t136 * t180 + t176) * t232 + (t241 * t172 - t203) * t151;
t239 = t131 * t141 * t142;
t166 = t171 * t167;
t224 = t168 * t171;
t190 = t176 * (t166 * t168 * t169 + t224);
t225 = t167 * t180;
t195 = t213 * t225;
t238 = (t169 * t195 + t177 * t190) / t164 ^ 2;
t237 = t142 * t171;
t236 = t142 * t182;
t197 = -qJD(1) * t172 + qJD(6);
t198 = qJD(6) * t172 - qJD(1);
t220 = t176 * t182;
t204 = t171 * t220;
t219 = t179 * t182;
t146 = -t198 * t219 + (t197 * t180 - t204) * t181;
t235 = t146 * t155 * t156;
t191 = t172 * t217 + t215;
t145 = t191 * qJD(1) - t163 * qJD(6) + t179 * t204;
t216 = t180 * t181;
t162 = t172 * t219 - t216;
t154 = t162 ^ 2;
t150 = t154 * t156 + 0.1e1;
t230 = t156 * t162;
t234 = 0.1e1 / t150 ^ 2 * (-t145 * t230 - t154 * t235);
t233 = t151 * t180;
t231 = t155 * t179;
t229 = t162 * t181;
t228 = t167 * t168;
t178 = t182 ^ 2;
t226 = t167 * t178;
t223 = t171 * t182;
t222 = t172 * t176;
t214 = qJD(1) * t180;
t139 = t142 * t226 + 0.1e1;
t212 = 0.2e1 * (-t226 * t239 + (t171 * t178 * t222 - t195) * t142) / t139 ^ 2;
t211 = 0.2e1 * t239;
t210 = -0.2e1 * t234;
t209 = t142 * t223;
t208 = t162 * t235;
t202 = 0.1e1 + t227;
t201 = t171 * t212;
t200 = -0.2e1 * t171 * t238;
t199 = t238 * t242;
t196 = t152 * t158 * t228;
t194 = t202 * t182;
t193 = t197 * t182;
t192 = t156 * t229 - t231;
t161 = -t172 * t216 + t219;
t148 = 0.1e1 / t150;
t147 = t202 * t180 * t158;
t137 = 0.1e1 / t139;
t135 = (t240 * t171 * t151 - t180 * t196) * t182;
t133 = -t172 * t233 + t232 + (t151 * t172 - t152 * t218) * t147;
t132 = -t202 * t199 + (qJD(1) * t194 + t190 * t242) * t158;
t129 = t192 * t210 * t223 + (t192 * t172 * t220 + (-t192 * t214 + ((-qJD(6) * t155 - 0.2e1 * t208) * t181 + (-t145 * t181 + (-qJD(6) * t162 + t146) * t179) * t156) * t182) * t171) * t148;
t128 = (t133 * t237 - t141 * t172) * t182 * t212 + ((-t141 * t214 + (-t133 * t176 - t131) * t236) * t172 + (-t141 * t220 - (-t132 * t152 * t180 - t241 * t151 + (t136 * t233 - t151 * t176 - t152 * t213) * t147) * t209 + (t142 * t214 + t182 * t211) * t133 - ((t132 - t213) * t151 + ((-t147 * t180 + 0.1e1) * t176 + (t147 - t180) * t136) * t152) * t172 * t236) * t171) * t137;
t1 = [t182 * t168 * t200 + (t176 * t194 - t214 * t224) * t158, 0, t132, t132, t132, 0; (t141 * t201 + (-t141 * t222 + (qJD(1) * t135 + t131) * t237) * t137) * t180 + (t142 * t201 * t135 + (-((t200 - t222 + (t136 * t168 * t225 + t222) * t158) * t151 + (t199 * t228 - t136 * t171 + (-t166 * t206 + (t136 - 0.2e1 * t221) * t171) * t158) * t152) * t209 + (-t142 * t222 + t171 * t211) * t135 + (-t141 + ((-t177 + t178) * t196 + t240 * t207) * t142) * t171 * qJD(1)) * t137) * t182, 0, t128, t128, t128, 0; 0.2e1 * (t155 * t191 + t161 * t230) * t234 + (0.2e1 * t161 * t208 - t198 * t155 * t216 + (t176 * t218 + t193) * t231 + (t161 * t145 + t191 * t146 - t193 * t229 - (t171 * t176 * t181 + t198 * t179) * t162 * t180) * t156) * t148, 0, t129, t129, t129, t210 + 0.2e1 * (-t145 * t156 * t148 + (-t148 * t235 - t156 * t234) * t162) * t162;];
JaD_rot  = t1;
