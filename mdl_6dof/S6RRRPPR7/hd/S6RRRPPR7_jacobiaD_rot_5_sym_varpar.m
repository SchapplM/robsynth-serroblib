% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:11
% EndTime: 2019-02-26 22:07:12
% DurationCPUTime: 0.81s
% Computational Cost: add. (975->110), mult. (3480->241), div. (475->12), fcn. (4141->11), ass. (0->108)
t169 = sin(qJ(2));
t160 = t169 ^ 2;
t172 = cos(qJ(2));
t163 = 0.1e1 / t172 ^ 2;
t225 = t160 * t163;
t170 = sin(qJ(1));
t161 = t170 ^ 2;
t155 = t161 * t225 + 0.1e1;
t162 = 0.1e1 / t172;
t222 = t162 * t169;
t241 = t169 * t225;
t181 = qJD(2) * (t162 * t241 + t222);
t173 = cos(qJ(1));
t213 = qJD(1) * t173;
t223 = t160 * t170;
t188 = t213 * t223;
t228 = (t161 * t181 + t163 * t188) / t155 ^ 2;
t242 = -0.2e1 * t228;
t196 = 0.1e1 + t225;
t240 = t170 * t196;
t168 = sin(qJ(3));
t171 = cos(qJ(3));
t209 = qJD(3) * t173;
t197 = t171 * t209;
t214 = qJD(1) * t170;
t239 = -t168 * t214 + t197;
t238 = t168 * t209 + t171 * t214;
t215 = t172 * t173;
t217 = t170 * t171;
t148 = t168 * t215 - t217;
t219 = t170 * t168;
t149 = t171 * t215 + t219;
t166 = sin(pkin(10));
t167 = cos(pkin(10));
t136 = t148 * t166 + t149 * t167;
t128 = 0.1e1 / t136;
t218 = t170 * t169;
t154 = atan2(t218, t172);
t151 = cos(t154);
t150 = sin(t154);
t204 = t150 * t218;
t140 = t151 * t172 + t204;
t137 = 0.1e1 / t140;
t129 = 0.1e1 / t136 ^ 2;
t138 = 0.1e1 / t140 ^ 2;
t237 = 0.2e1 * t169;
t152 = 0.1e1 / t155;
t236 = t152 - 0.1e1;
t165 = t173 ^ 2;
t224 = t160 * t165;
t126 = t138 * t224 + 0.1e1;
t211 = qJD(2) * t172;
t202 = t169 * t213;
t212 = qJD(2) * t170;
t118 = ((t170 * t211 + t202) * t162 + t212 * t225) * t152;
t226 = t151 * t169;
t109 = (t118 * t170 - qJD(2)) * t226 + (t202 + (-t118 + t212) * t172) * t150;
t234 = t109 * t137 * t138;
t235 = (-t224 * t234 + (t165 * t169 * t211 - t188) * t138) / t126 ^ 2;
t233 = t118 * t150;
t232 = t118 * t169;
t192 = -t148 * t167 + t149 * t166;
t231 = t129 * t192;
t184 = -t166 * t168 - t167 * t171;
t220 = t169 * t173;
t144 = t184 * t220;
t230 = t129 * t144;
t229 = t138 * t169;
t142 = t152 * t240;
t227 = t142 * t170;
t221 = t168 * t173;
t216 = t170 * t172;
t210 = qJD(2) * t173;
t182 = t168 * t216 + t171 * t173;
t199 = t169 * t210;
t120 = t182 * qJD(1) - qJD(3) * t219 + t168 * t199 - t172 * t197;
t190 = -qJD(1) * t172 + qJD(3);
t191 = qJD(3) * t172 - qJD(1);
t121 = -t191 * t221 + (t190 * t170 - t199) * t171;
t112 = t120 * t167 + t121 * t166;
t127 = t192 ^ 2;
t116 = t127 * t129 + 0.1e1;
t130 = t128 * t129;
t186 = t120 * t166 - t121 * t167;
t208 = 0.2e1 * (t127 * t130 * t186 + t112 * t231) / t116 ^ 2;
t207 = -0.2e1 * t234;
t206 = 0.2e1 * t130 * t192;
t205 = t138 * t220;
t203 = t152 * t160 * t162;
t195 = -0.2e1 * t169 * t235;
t194 = t186 * t206;
t193 = t162 * t242;
t189 = t170 * t203;
t187 = t196 * t173;
t185 = -t166 * t171 + t167 * t168;
t183 = t190 * t173;
t147 = -t171 * t216 + t221;
t143 = t185 * t220;
t132 = t147 * t167 - t166 * t182;
t131 = t147 * t166 + t167 * t182;
t124 = 0.1e1 / t126;
t123 = t171 * t183 + (qJD(2) * t169 * t171 + t191 * t168) * t170;
t122 = -t191 * t217 + (t169 * t212 + t183) * t168;
t117 = (-t236 * t169 * t150 + t151 * t189) * t173;
t114 = 0.1e1 / t116;
t111 = t150 * t216 - t226 + (-t150 * t172 + t151 * t218) * t142;
t110 = t240 * t242 + (qJD(1) * t187 + 0.2e1 * t170 * t181) * t152;
t1 = [t193 * t220 + (qJD(2) * t187 - t214 * t222) * t152, t110, 0, 0, 0, 0; (t137 * t195 + (t137 * t211 + (-qJD(1) * t117 - t109) * t229) * t124) * t170 + (t138 * t195 * t117 + (((-t118 * t189 - t236 * t211 + t228 * t237) * t150 + (t193 * t223 + t232 + (-t232 + (t237 + t241) * t212) * t152) * t151) * t205 + (t138 * t211 + t169 * t207) * t117 + (t137 + ((-t161 + t165) * t151 * t203 + t236 * t204) * t138) * t169 * qJD(1)) * t124) * t173, 0.2e1 * (-t111 * t229 + t137 * t172) * t173 * t235 + ((t137 * t214 + (qJD(2) * t111 + t109) * t173 * t138) * t172 + (t137 * t210 + (t110 * t151 * t170 - t150 * t212 - t227 * t233 + t233 + (qJD(2) * t150 + t151 * t213) * t142) * t205 + (-t138 * t214 + t173 * t207) * t111 + ((-t110 + t213) * t150 + ((-0.1e1 + t227) * qJD(2) + (-t142 + t170) * t118) * t151) * t138 * t215) * t169) * t124, 0, 0, 0, 0; (-t128 * t131 + t132 * t231) * t208 + ((-t122 * t167 + t123 * t166) * t128 - t132 * t194 + (t131 * t186 - (t122 * t166 + t123 * t167) * t192 - t132 * t112) * t129) * t114 (-t128 * t143 + t192 * t230) * t208 + (-t112 * t230 - (-t129 * t143 + t144 * t206) * t186 + (t185 * t128 - t184 * t231) * t172 * t210 + ((t239 * t128 - t238 * t231) * t167 + (t238 * t128 + t239 * t231) * t166) * t169) * t114 (t128 * t136 + t192 * t231) * t208 + (t186 * t128 - t192 * t194 + (-0.2e1 * t192 * t112 - t136 * t186) * t129) * t114, 0, 0, 0;];
JaD_rot  = t1;
