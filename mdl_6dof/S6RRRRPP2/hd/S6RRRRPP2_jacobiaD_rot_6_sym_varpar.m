% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP2
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
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:25:31
% EndTime: 2019-02-26 22:25:32
% DurationCPUTime: 0.75s
% Computational Cost: add. (3091->99), mult. (3836->212), div. (757->12), fcn. (4483->9), ass. (0->95)
t173 = sin(qJ(1));
t169 = t173 ^ 2;
t171 = qJ(2) + qJ(3);
t166 = sin(t171);
t162 = t166 ^ 2;
t167 = cos(t171);
t164 = 0.1e1 / t167 ^ 2;
t216 = t162 * t164;
t158 = t169 * t216 + 0.1e1;
t161 = t166 * t162;
t163 = 0.1e1 / t167;
t168 = qJD(2) + qJD(3);
t213 = t163 * t166;
t183 = t168 * (t161 * t163 * t164 + t213);
t175 = cos(qJ(1));
t202 = qJD(1) * t175;
t214 = t162 * t173;
t186 = t202 * t214;
t225 = (t164 * t186 + t169 * t183) / t158 ^ 2;
t232 = -0.2e1 * t225;
t172 = sin(qJ(4));
t174 = cos(qJ(4));
t209 = t168 * t175;
t194 = t166 * t209;
t203 = qJD(1) * t173;
t205 = t173 * t174;
t208 = t172 * t175;
t135 = t174 * t194 - qJD(4) * t205 - t172 * t202 + (qJD(4) * t208 + t174 * t203) * t167;
t151 = -t167 * t208 + t205;
t145 = t151 ^ 2;
t204 = t174 * t175;
t206 = t173 * t172;
t152 = t167 * t204 + t206;
t147 = 0.1e1 / t152 ^ 2;
t222 = t145 * t147;
t144 = 0.1e1 + t222;
t188 = qJD(1) * t167 - qJD(4);
t189 = qJD(4) * t167 - qJD(1);
t134 = -t189 * t204 + (t188 * t173 + t194) * t172;
t219 = t147 * t151;
t199 = t134 * t219;
t146 = 0.1e1 / t152;
t148 = t146 * t147;
t221 = t145 * t148;
t231 = (t135 * t221 + t199) / t144 ^ 2;
t192 = 0.1e1 + t216;
t230 = t173 * t192;
t155 = 0.1e1 / t158;
t193 = t166 * t202;
t210 = t168 * t173;
t195 = t164 * t210;
t130 = ((t167 * t210 + t193) * t163 + t162 * t195) * t155;
t229 = t130 - t210;
t207 = t173 * t166;
t157 = atan2(t207, t167);
t154 = cos(t157);
t153 = sin(t157);
t196 = t153 * t207;
t140 = t154 * t167 + t196;
t137 = 0.1e1 / t140;
t138 = 0.1e1 / t140 ^ 2;
t228 = t155 - 0.1e1;
t170 = t175 ^ 2;
t215 = t162 * t170;
t133 = t138 * t215 + 0.1e1;
t211 = t167 * t168;
t217 = t154 * t166;
t125 = (t130 * t173 - t168) * t217 + (-t167 * t229 + t193) * t153;
t226 = t125 * t137 * t138;
t227 = (-t215 * t226 + (t166 * t170 * t211 - t186) * t138) / t133 ^ 2;
t224 = t138 * t166;
t223 = t138 * t175;
t220 = t146 * t172;
t218 = t153 * t173;
t212 = t166 * t175;
t201 = -0.2e1 * t226;
t200 = 0.2e1 * t231;
t198 = t135 * t148 * t151;
t197 = t138 * t212;
t191 = -0.2e1 * t166 * t227;
t190 = t163 * t232;
t187 = t154 * t155 * t162 * t163;
t185 = t192 * t175;
t184 = t174 * t219 + t220;
t150 = -t167 * t205 + t208;
t149 = t167 * t206 + t204;
t142 = 0.1e1 / t144;
t141 = t155 * t230;
t131 = 0.1e1 / t133;
t129 = (-t228 * t166 * t153 + t173 * t187) * t175;
t127 = t167 * t218 - t217 + (-t153 * t167 + t154 * t207) * t141;
t126 = t230 * t232 + (qJD(1) * t185 + 0.2e1 * t173 * t183) * t155;
t123 = -0.2e1 * t184 * t212 * t231 + (t184 * t167 * t209 + (-t184 * t203 + ((qJD(4) * t146 + 0.2e1 * t198) * t174 + (t134 * t174 + (-qJD(4) * t151 + t135) * t172) * t147) * t175) * t166) * t142;
t122 = 0.2e1 * (-t127 * t224 + t137 * t167) * t175 * t227 + ((t137 * t203 + (t127 * t168 + t125) * t223) * t167 + (t137 * t209 + (t126 * t154 * t173 + t229 * t153 + (-t130 * t218 + t153 * t168 + t154 * t202) * t141) * t197 + (-t138 * t203 + t175 * t201) * t127 + ((-t126 + t202) * t153 + ((t141 * t173 - 0.1e1) * t168 + (-t141 + t173) * t130) * t154) * t167 * t223) * t166) * t131;
t1 = [t190 * t212 + (t168 * t185 - t203 * t213) * t155, t126, t126, 0, 0, 0; (t137 * t191 + (t137 * t211 + (-qJD(1) * t129 - t125) * t224) * t131) * t173 + (t138 * t191 * t129 + (((0.2e1 * t166 * t225 + t211 + (-t130 * t163 * t214 - t211) * t155) * t153 + (t190 * t214 + t130 * t166 + (t161 * t195 + (-t130 + 0.2e1 * t210) * t166) * t155) * t154) * t197 + (t138 * t211 + t166 * t201) * t129 + (t137 + ((-t169 + t170) * t187 + t228 * t196) * t138) * t166 * qJD(1)) * t131) * t175, t122, t122, 0, 0, 0; (-t146 * t149 + t150 * t219) * t200 + (-0.2e1 * t150 * t198 + t189 * t146 * t205 + (-t168 * t207 + t188 * t175) * t220 + (-t150 * t134 + t149 * t135 + (t188 * t204 - (t166 * t168 * t174 + t189 * t172) * t173) * t151) * t147) * t142, t123, t123 (t146 * t152 + t222) * t200 + (-0.2e1 * t199 + (-t147 * t152 + t146 - 0.2e1 * t221) * t135) * t142, 0, 0;];
JaD_rot  = t1;
