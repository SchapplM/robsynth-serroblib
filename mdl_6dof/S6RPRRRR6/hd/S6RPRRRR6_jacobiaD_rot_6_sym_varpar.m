% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR6
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
% Datum: 2019-02-26 21:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:17:36
% EndTime: 2019-02-26 21:17:37
% DurationCPUTime: 0.78s
% Computational Cost: add. (3714->96), mult. (2949->203), div. (516->12), fcn. (3430->9), ass. (0->95)
t169 = pkin(11) + qJ(3);
t165 = sin(t169);
t161 = t165 ^ 2;
t166 = cos(t169);
t163 = 0.1e1 / t166 ^ 2;
t216 = t161 * t163;
t172 = sin(qJ(1));
t234 = 0.2e1 * t172;
t233 = t165 * t216;
t170 = t172 ^ 2;
t155 = t170 * t216 + 0.1e1;
t153 = 0.1e1 / t155;
t162 = 0.1e1 / t166;
t173 = cos(qJ(1));
t207 = qJD(1) * t173;
t195 = t165 * t207;
t205 = qJD(3) * t172;
t128 = (-(-t166 * t205 - t195) * t162 + t205 * t216) * t153;
t232 = t128 - t205;
t168 = qJ(4) + qJ(5) + qJ(6);
t159 = cos(t168);
t158 = sin(t168);
t211 = t172 * t158;
t212 = t166 * t173;
t148 = t159 * t212 + t211;
t209 = t172 * t165;
t152 = atan2(-t209, -t166);
t151 = cos(t152);
t150 = sin(t152);
t198 = t150 * t209;
t137 = -t151 * t166 - t198;
t132 = 0.1e1 / t137;
t142 = 0.1e1 / t148;
t133 = 0.1e1 / t137 ^ 2;
t143 = 0.1e1 / t148 ^ 2;
t231 = -0.2e1 * t165;
t230 = t153 - 0.1e1;
t219 = t151 * t165;
t121 = (-t128 * t172 + qJD(3)) * t219 + (t166 * t232 - t195) * t150;
t229 = t121 * t132 * t133;
t167 = qJD(4) + qJD(5) + qJD(6);
t183 = t159 * t173 + t166 * t211;
t204 = qJD(3) * t173;
t194 = t165 * t204;
t126 = t183 * qJD(1) - t148 * t167 + t158 * t194;
t210 = t172 * t159;
t147 = t158 * t212 - t210;
t141 = t147 ^ 2;
t138 = t141 * t143 + 0.1e1;
t222 = t143 * t147;
t188 = -qJD(1) * t166 + t167;
t189 = t166 * t167 - qJD(1);
t218 = t158 * t173;
t127 = -t189 * t218 + (t188 * t172 - t194) * t159;
t227 = t127 * t142 * t143;
t228 = (-t126 * t222 - t141 * t227) / t138 ^ 2;
t226 = t128 * t165;
t225 = t133 * t165;
t214 = t162 * t165;
t182 = qJD(3) * (t162 * t233 + t214);
t186 = t161 * t172 * t207;
t224 = (t163 * t186 + t170 * t182) / t155 ^ 2;
t223 = t142 * t158;
t221 = t147 * t159;
t220 = t150 * t172;
t217 = t161 * t162;
t171 = t173 ^ 2;
t215 = t161 * t171;
t213 = t165 * t173;
t208 = qJD(1) * t172;
t206 = qJD(3) * t166;
t131 = t133 * t215 + 0.1e1;
t203 = 0.2e1 * (-t215 * t229 + (t165 * t171 * t206 - t186) * t133) / t131 ^ 2;
t202 = 0.2e1 * t229;
t201 = -0.2e1 * t228;
t200 = t147 * t227;
t199 = t133 * t213;
t197 = t153 * t217;
t193 = 0.1e1 + t216;
t192 = t165 * t203;
t191 = t224 * t231;
t190 = t224 * t234;
t187 = t172 * t197;
t185 = t193 * t173;
t184 = t143 * t221 - t223;
t181 = t165 * t205 + t188 * t173;
t146 = -t166 * t210 + t218;
t140 = t193 * t172 * t153;
t135 = 0.1e1 / t138;
t129 = 0.1e1 / t131;
t125 = (t230 * t165 * t150 - t151 * t187) * t173;
t124 = -t166 * t220 + t219 + (t150 * t166 - t151 * t209) * t140;
t122 = -t193 * t190 + (qJD(1) * t185 + t182 * t234) * t153;
t119 = t201 + 0.2e1 * (-t126 * t135 * t143 + (-t135 * t227 - t143 * t228) * t147) * t147;
t1 = [t162 * t173 * t191 + (qJD(3) * t185 - t208 * t214) * t153, 0, t122, 0, 0, 0; (t132 * t192 + (-t132 * t206 + (qJD(1) * t125 + t121) * t225) * t129) * t172 + (t133 * t192 * t125 + (-((t128 * t187 + t230 * t206 + t191) * t150 + (t190 * t217 - t226 + (t226 + (t231 - t233) * t205) * t153) * t151) * t199 + (-t133 * t206 + t165 * t202) * t125 + (-t132 + ((-t170 + t171) * t151 * t197 + t230 * t198) * t133) * t165 * qJD(1)) * t129) * t173, 0 (t124 * t225 - t132 * t166) * t173 * t203 + ((-t132 * t208 + (-qJD(3) * t124 - t121) * t173 * t133) * t166 + (-t132 * t204 - (-t122 * t151 * t172 - t232 * t150 + (-qJD(3) * t150 + t128 * t220 - t151 * t207) * t140) * t199 + (t133 * t208 + t173 * t202) * t124 - ((t122 - t207) * t150 + ((-t140 * t172 + 0.1e1) * qJD(3) + (t140 - t172) * t128) * t151) * t133 * t212) * t165) * t129, 0, 0, 0; 0.2e1 * (t142 * t183 + t146 * t222) * t228 + (0.2e1 * t146 * t200 - t189 * t142 * t210 + t181 * t223 + (-t189 * t147 * t211 + t146 * t126 + t127 * t183 - t181 * t221) * t143) * t135, 0, t184 * t201 * t213 + (t184 * t166 * t204 + (-t184 * t208 + ((-t142 * t167 - 0.2e1 * t200) * t159 + (-t126 * t159 + (-t147 * t167 + t127) * t158) * t143) * t173) * t165) * t135, t119, t119, t119;];
JaD_rot  = t1;
