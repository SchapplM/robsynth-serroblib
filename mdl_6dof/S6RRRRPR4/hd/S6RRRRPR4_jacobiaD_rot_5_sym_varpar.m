% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:32:27
% EndTime: 2019-02-26 22:32:28
% DurationCPUTime: 0.73s
% Computational Cost: add. (4131->97), mult. (3810->206), div. (753->12), fcn. (4455->9), ass. (0->95)
t172 = sin(qJ(1));
t232 = 0.2e1 * t172;
t169 = t172 ^ 2;
t171 = qJ(2) + qJ(3);
t165 = sin(t171);
t159 = t165 ^ 2;
t166 = cos(t171);
t161 = 0.1e1 / t166 ^ 2;
t218 = t159 * t161;
t155 = t169 * t218 + 0.1e1;
t152 = 0.1e1 / t155;
t160 = 0.1e1 / t166;
t173 = cos(qJ(1));
t204 = qJD(1) * t173;
t194 = t165 * t204;
t167 = qJD(2) + qJD(3);
t210 = t167 * t172;
t197 = t161 * t210;
t126 = (-(-t166 * t210 - t194) * t160 + t159 * t197) * t152;
t231 = t126 - t210;
t168 = qJ(4) + pkin(11);
t164 = cos(t168);
t163 = sin(t168);
t208 = t172 * t163;
t211 = t166 * t173;
t148 = t164 * t211 + t208;
t206 = t172 * t165;
t151 = atan2(-t206, -t166);
t150 = cos(t151);
t149 = sin(t151);
t198 = t149 * t206;
t136 = -t150 * t166 - t198;
t133 = 0.1e1 / t136;
t142 = 0.1e1 / t148;
t134 = 0.1e1 / t136 ^ 2;
t143 = 0.1e1 / t148 ^ 2;
t230 = t152 - 0.1e1;
t220 = t150 * t165;
t121 = (-t126 * t172 + t167) * t220 + (t231 * t166 - t194) * t149;
t229 = t121 * t133 * t134;
t183 = t164 * t173 + t166 * t208;
t209 = t167 * t173;
t195 = t165 * t209;
t127 = t183 * qJD(1) - t148 * qJD(4) + t163 * t195;
t207 = t172 * t164;
t147 = t163 * t211 - t207;
t141 = t147 ^ 2;
t139 = t141 * t143 + 0.1e1;
t223 = t143 * t147;
t188 = -qJD(1) * t166 + qJD(4);
t189 = qJD(4) * t166 - qJD(1);
t214 = t163 * t173;
t128 = -t189 * t214 + (t172 * t188 - t195) * t164;
t227 = t128 * t142 * t143;
t228 = (-t127 * t223 - t141 * t227) / t139 ^ 2;
t158 = t165 * t159;
t215 = t160 * t165;
t182 = t167 * (t158 * t160 * t161 + t215);
t216 = t159 * t172;
t186 = t204 * t216;
t226 = (t161 * t186 + t169 * t182) / t155 ^ 2;
t225 = t134 * t165;
t224 = t142 * t163;
t222 = t147 * t164;
t221 = t149 * t172;
t219 = t159 * t160;
t170 = t173 ^ 2;
t217 = t159 * t170;
t213 = t165 * t173;
t212 = t166 * t167;
t205 = qJD(1) * t172;
t131 = t134 * t217 + 0.1e1;
t203 = 0.2e1 * (-t217 * t229 + (t165 * t170 * t212 - t186) * t134) / t131 ^ 2;
t202 = 0.2e1 * t229;
t201 = -0.2e1 * t228;
t200 = t134 * t213;
t199 = t147 * t227;
t193 = 0.1e1 + t218;
t192 = t165 * t203;
t191 = -0.2e1 * t165 * t226;
t190 = t226 * t232;
t187 = t150 * t152 * t219;
t185 = t193 * t173;
t184 = t143 * t222 - t224;
t181 = t167 * t206 + t173 * t188;
t146 = -t166 * t207 + t214;
t140 = t193 * t172 * t152;
t137 = 0.1e1 / t139;
t129 = 0.1e1 / t131;
t125 = (t149 * t165 * t230 - t172 * t187) * t173;
t124 = -t166 * t221 + t220 + (t149 * t166 - t150 * t206) * t140;
t122 = -t193 * t190 + (qJD(1) * t185 + t182 * t232) * t152;
t119 = t184 * t201 * t213 + (t184 * t166 * t209 + (-t184 * t205 + ((-qJD(4) * t142 - 0.2e1 * t199) * t164 + (-t127 * t164 + (-qJD(4) * t147 + t128) * t163) * t143) * t173) * t165) * t137;
t118 = (t124 * t225 - t133 * t166) * t173 * t203 + ((-t133 * t205 + (-t124 * t167 - t121) * t173 * t134) * t166 + (-t133 * t209 - (-t122 * t150 * t172 - t231 * t149 + (t126 * t221 - t149 * t167 - t150 * t204) * t140) * t200 + (t134 * t205 + t173 * t202) * t124 - ((t122 - t204) * t149 + ((-t140 * t172 + 0.1e1) * t167 + (t140 - t172) * t126) * t150) * t134 * t211) * t165) * t129;
t1 = [t160 * t173 * t191 + (t167 * t185 - t205 * t215) * t152, t122, t122, 0, 0, 0; (t133 * t192 + (-t133 * t212 + (qJD(1) * t125 + t121) * t225) * t129) * t172 + (t134 * t192 * t125 + (-((t191 - t212 + (t126 * t160 * t216 + t212) * t152) * t149 + (t190 * t219 - t126 * t165 + (-t158 * t197 + (t126 - 0.2e1 * t210) * t165) * t152) * t150) * t200 + (-t134 * t212 + t165 * t202) * t125 + (-t133 + ((-t169 + t170) * t187 + t230 * t198) * t134) * t165 * qJD(1)) * t129) * t173, t118, t118, 0, 0, 0; 0.2e1 * (t142 * t183 + t146 * t223) * t228 + (0.2e1 * t146 * t199 - t189 * t142 * t207 + t181 * t224 + (-t147 * t189 * t208 + t146 * t127 + t128 * t183 - t181 * t222) * t143) * t137, t119, t119, t201 + 0.2e1 * (-t127 * t143 * t137 + (-t137 * t227 - t143 * t228) * t147) * t147, 0, 0;];
JaD_rot  = t1;
