% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:24
% EndTime: 2019-02-26 22:28:25
% DurationCPUTime: 1.04s
% Computational Cost: add. (6849->120), mult. (8141->260), div. (1608->14), fcn. (10341->9), ass. (0->113)
t188 = qJ(3) + qJ(4);
t178 = sin(t188);
t179 = cos(t188);
t192 = cos(qJ(1));
t237 = t192 * t179;
t190 = sin(qJ(1));
t191 = cos(qJ(2));
t239 = t190 * t191;
t160 = t178 * t239 + t237;
t180 = qJD(3) + qJD(4);
t266 = t160 * t180;
t189 = sin(qJ(2));
t265 = t189 * t192;
t164 = t190 * t178 + t191 * t237;
t181 = 0.1e1 / t189;
t182 = 0.1e1 / t189 ^ 2;
t183 = t181 * t182;
t264 = qJD(2) * (0.2e1 * t183 * t191 ^ 2 + t181);
t232 = qJD(2) * t192;
t212 = t189 * t232;
t144 = t160 * qJD(1) - t164 * t180 + t178 * t212;
t238 = t192 * t178;
t240 = t190 * t179;
t163 = t191 * t238 - t240;
t175 = 0.1e1 / t178;
t176 = 0.1e1 / t178 ^ 2;
t233 = qJD(2) * t191;
t215 = t182 * t233;
t248 = t179 * t180;
t250 = t175 * t181;
t263 = (t176 * t181 * t248 + t175 * t215) * t163 + t144 * t250;
t243 = t189 * t178;
t152 = atan2(-t160, t243);
t149 = cos(t152);
t148 = sin(t152);
t256 = t148 * t160;
t143 = t149 * t243 - t256;
t140 = 0.1e1 / t143;
t185 = 0.1e1 / t192;
t141 = 0.1e1 / t143 ^ 2;
t186 = 0.1e1 / t192 ^ 2;
t262 = -0.2e1 * t160;
t261 = 0.2e1 * t163;
t157 = t160 ^ 2;
t249 = t176 * t182;
t153 = t157 * t249 + 0.1e1;
t150 = 0.1e1 / t153;
t247 = t179 * t189;
t203 = t178 * t233 + t180 * t247;
t222 = t160 * t249;
t242 = t189 * t190;
t213 = qJD(2) * t242;
t234 = qJD(1) * t192;
t235 = qJD(1) * t190;
t146 = -t178 * t213 - t180 * t238 - t179 * t235 + (t178 * t234 + t180 * t240) * t191;
t224 = t146 * t250;
t132 = (t203 * t222 - t224) * t150;
t201 = -t132 * t160 + t203;
t127 = (-t132 * t243 - t146) * t148 + t201 * t149;
t142 = t140 * t141;
t260 = t127 * t142;
t177 = t175 * t176;
t214 = t183 * t233;
t220 = t182 * t248;
t259 = (t146 * t222 + (-t176 * t214 - t177 * t220) * t157) / t153 ^ 2;
t258 = t141 * t163;
t257 = t144 * t141;
t255 = t148 * t163;
t254 = t148 * t189;
t253 = t149 * t160;
t252 = t149 * t163;
t251 = t149 * t191;
t246 = t182 * t186;
t245 = t182 * t191;
t244 = t186 * t190;
t221 = t175 * t245;
t206 = t160 * t221 + t190;
t139 = t206 * t150;
t236 = -t139 + t190;
t158 = t163 ^ 2;
t138 = t158 * t141 + 0.1e1;
t231 = 0.2e1 * (-t158 * t260 - t163 * t257) / t138 ^ 2;
t230 = -0.2e1 * t259;
t145 = (-t180 * t191 + qJD(1)) * t238 + (-t212 + (-qJD(1) * t191 + t180) * t190) * t179;
t159 = t164 ^ 2;
t156 = t159 * t246 + 0.1e1;
t187 = t185 * t186;
t229 = 0.2e1 * (t164 * t145 * t246 + (t182 * t187 * t235 - t186 * t214) * t159) / t156 ^ 2;
t228 = t142 * t261;
t227 = t181 * t259;
t226 = t141 * t255;
t223 = t160 * t250;
t219 = t185 * t245;
t216 = t186 * t235;
t211 = t140 * t231;
t210 = t141 * t231;
t209 = t181 * t229;
t207 = t175 * t227;
t162 = t179 * t239 - t238;
t205 = t160 * t176 * t179 - t162 * t175;
t204 = t162 * t185 - t164 * t244;
t154 = 0.1e1 / t156;
t147 = t164 * qJD(1) - t179 * t213 - t266;
t136 = 0.1e1 / t138;
t135 = t205 * t181 * t150;
t131 = (-t148 + (t149 * t223 + t148) * t150) * t163;
t130 = -t139 * t253 + (t236 * t254 + t251) * t178;
t129 = t163 * t185 * t209 + (t144 * t181 * t185 + (-t181 * t216 + t185 * t215) * t163) * t154;
t128 = t149 * t247 - t148 * t162 + (-t148 * t243 - t253) * t135;
t126 = t206 * t230 + (t146 * t221 + t234 + (-t176 * t191 * t220 - t175 * t264) * t160) * t150;
t124 = -0.2e1 * t205 * t227 + (-t205 * t215 + ((-t147 - t266) * t175 + (t177 * t248 * t262 + (t162 * t180 + t146) * t176) * t179) * t181) * t150;
t123 = (t128 * t258 - t140 * t164) * t231 + (t128 * t257 + t145 * t140 + (t128 * t228 - t164 * t141) * t127 - (t179 * t233 - t180 * t243 - t124 * t160 - t135 * t146 + (-t135 * t243 - t162) * t132) * t141 * t252 - (-t147 + (-t124 * t178 - t132 * t179) * t189 - t201 * t135) * t226) * t136;
t1 = [t263 * t150 + t207 * t261, t126, t124, t124, 0, 0; t160 * t211 + (-t146 * t140 + (t127 * t160 + t131 * t144) * t141) * t136 + (t131 * t210 + (0.2e1 * t131 * t260 + (t144 * t150 - t144 - (-t132 * t150 * t223 + t230) * t163) * t141 * t148 + (-(t207 * t262 - t132) * t258 + (-(t132 + t224) * t163 + t263 * t160) * t141 * t150) * t149) * t136) * t163, t130 * t163 * t210 + (-(-t126 * t253 + (t132 * t256 - t146 * t149) * t139) * t258 + (-t140 * t265 - (-t139 * t254 + t148 * t242 + t251) * t258) * t248 + (t127 * t228 + t257) * t130) * t136 + (t211 * t265 + ((-t140 * t232 - (t236 * qJD(2) - t132) * t226) * t191 + (t140 * t235 + (t192 * t127 - (-t126 + t234) * t255 - (t236 * t132 - qJD(2)) * t252) * t141) * t189) * t136) * t178, t123, t123, 0, 0; t204 * t209 + (t204 * t215 + (t145 * t244 - t147 * t185 + (-t162 * t244 + (0.2e1 * t187 * t190 ^ 2 + t185) * t164) * qJD(1)) * t181) * t154 (t164 * t219 + t179) * t229 + (-t145 * t219 + t178 * t180 + (t185 * t264 - t216 * t245) * t164) * t154, t129, t129, 0, 0;];
JaD_rot  = t1;
