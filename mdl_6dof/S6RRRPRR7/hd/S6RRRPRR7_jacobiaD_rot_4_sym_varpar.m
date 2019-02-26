% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR7
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR7_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobiaD_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:18
% EndTime: 2019-02-26 22:19:19
% DurationCPUTime: 0.72s
% Computational Cost: add. (1788->92), mult. (4303->201), div. (668->14), fcn. (5516->11), ass. (0->91)
t179 = sin(qJ(2));
t180 = sin(qJ(1));
t231 = cos(pkin(6));
t202 = t180 * t231;
t200 = t179 * t202;
t181 = cos(qJ(2));
t182 = cos(qJ(1));
t216 = t182 * t181;
t163 = -t200 + t216;
t174 = qJ(3) + pkin(12);
t170 = sin(t174);
t171 = cos(t174);
t178 = sin(pkin(6));
t220 = t178 * t180;
t192 = -t163 * t170 + t171 * t220;
t233 = t192 * qJD(3);
t201 = t182 * t231;
t199 = t181 * t201;
t217 = t180 * t179;
t159 = -t199 + t217;
t219 = t178 * t181;
t153 = atan2(-t159, -t219);
t151 = sin(t153);
t152 = cos(t153);
t157 = t159 ^ 2;
t173 = 0.1e1 / t178 ^ 2;
t176 = 0.1e1 / t181 ^ 2;
t156 = t157 * t173 * t176 + 0.1e1;
t154 = 0.1e1 / t156;
t172 = 0.1e1 / t178;
t175 = 0.1e1 / t181;
t206 = t159 * t172 * t175;
t232 = (t152 * t206 - t151) * t154 + t151;
t135 = -t151 * t159 - t152 * t219;
t132 = 0.1e1 / t135;
t150 = t163 * t171 + t170 * t220;
t144 = 0.1e1 / t150;
t133 = 0.1e1 / t135 ^ 2;
t145 = 0.1e1 / t150 ^ 2;
t189 = -t179 * t201 - t180 * t181;
t190 = -t182 * t179 - t181 * t202;
t141 = -t190 * qJD(1) - t189 * qJD(2);
t214 = qJD(2) * t179;
t203 = t176 * t214;
t191 = t141 * t175 + t159 * t203;
t222 = t154 * t172;
t124 = t191 * t222;
t195 = t151 * t219 - t152 * t159;
t207 = t152 * t178 * t179;
t120 = qJD(2) * t207 + t195 * t124 - t151 * t141;
t230 = t120 * t132 * t133;
t140 = t189 * qJD(1) + t190 * qJD(2);
t215 = qJD(1) * t178;
t204 = t182 * t215;
t130 = t140 * t171 + t170 * t204 + t233;
t229 = t130 * t144 * t145;
t177 = t175 * t176;
t228 = (t141 * t159 * t176 + t157 * t177 * t214) * t173 / t156 ^ 2;
t227 = t133 * t190;
t129 = t150 * qJD(3) + t140 * t170 - t171 * t204;
t143 = t192 ^ 2;
t138 = t143 * t145 + 0.1e1;
t225 = t145 * t192;
t226 = 0.1e1 / t138 ^ 2 * (-t129 * t225 - t143 * t229);
t224 = t151 * t190;
t223 = t152 * t190;
t221 = t176 * t179;
t218 = t178 * t182;
t213 = qJD(2) * t181;
t158 = t190 ^ 2;
t128 = t133 * t158 + 0.1e1;
t198 = qJD(2) * t231 + qJD(1);
t139 = -qJD(1) * t199 - t182 * t213 + t198 * t217;
t212 = 0.2e1 * (t139 * t227 - t158 * t230) / t128 ^ 2;
t211 = 0.2e1 * t230;
t210 = -0.2e1 * t228;
t209 = 0.2e1 * t226;
t208 = t192 * t229;
t205 = t180 * t215;
t196 = t144 * t170 + t171 * t225;
t194 = t159 * t221 - t175 * t189;
t193 = -t170 * t189 + t171 * t218;
t148 = t170 * t218 + t171 * t189;
t142 = -qJD(1) * t200 - t180 * t214 + t198 * t216;
t136 = 0.1e1 / t138;
t126 = 0.1e1 / t128;
t125 = t194 * t222;
t123 = t232 * t190;
t121 = t195 * t125 + t151 * t189 + t207;
t119 = (t194 * t210 + (t141 * t221 + t142 * t175 + (-t189 * t221 + (0.2e1 * t177 * t179 ^ 2 + t175) * t159) * qJD(2)) * t154) * t172;
t1 = [(-t190 * t175 * t210 + (-t139 * t175 - t190 * t203) * t154) * t172, t119, 0, 0, 0, 0; t159 * t132 * t212 + (-t141 * t132 + (t120 * t159 + t123 * t139) * t133) * t126 - (t123 * t211 * t126 + (t123 * t212 + ((t124 * t154 * t206 + t210) * t224 + (0.2e1 * t206 * t228 - t124 + (-t191 * t172 + t124) * t154) * t223 - t232 * t139) * t126) * t133) * t190 (-t121 * t227 - t132 * t163) * t212 + (-t121 * t190 * t211 + t140 * t132 + (-t163 * t120 + t121 * t139 + (t178 * t213 - t119 * t159 - t125 * t141 + (t125 * t219 + t189) * t124) * t223 + (t124 * t125 * t159 - t142 + (t119 * t181 + (-qJD(2) * t125 - t124) * t179) * t178) * t224) * t133) * t126, 0, 0, 0, 0; (t144 * t193 - t148 * t225) * t209 + ((t148 * qJD(3) - t142 * t170 + t171 * t205) * t144 - 0.2e1 * t148 * t208 + (t193 * t130 + (t193 * qJD(3) - t142 * t171 - t170 * t205) * t192 - t148 * t129) * t145) * t136, -t196 * t190 * t209 + (t196 * t139 - ((-qJD(3) * t144 + 0.2e1 * t208) * t171 + (t129 * t171 + (t130 + t233) * t170) * t145) * t190) * t136, -0.2e1 * t226 - 0.2e1 * (t129 * t145 * t136 - (-t136 * t229 - t145 * t226) * t192) * t192, 0, 0, 0;];
JaD_rot  = t1;
