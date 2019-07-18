% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
%
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RRPRR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_jacobiaD_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_jacobiaD_rot_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_jacobiaD_rot_5_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:52
% EndTime: 2019-07-18 17:22:53
% DurationCPUTime: 0.75s
% Computational Cost: add. (3645->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
t172 = sin(qJ(1));
t234 = 0.2e1 * t172;
t168 = t172 ^ 2;
t170 = qJ(2) + qJ(4);
t165 = sin(t170);
t161 = t165 ^ 2;
t166 = cos(t170);
t163 = 0.1e1 / t166 ^ 2;
t219 = t161 * t163;
t156 = t168 * t219 + 0.1e1;
t154 = 0.1e1 / t156;
t162 = 0.1e1 / t166;
t174 = cos(qJ(1));
t205 = qJD(1) * t174;
t195 = t165 * t205;
t167 = qJD(2) + qJD(4);
t213 = t167 * t172;
t198 = t163 * t213;
t128 = (-(-t166 * t213 - t195) * t162 + t161 * t198) * t154;
t233 = t128 - t213;
t173 = cos(qJ(5));
t207 = t173 * t174;
t171 = sin(qJ(5));
t209 = t172 * t171;
t150 = t166 * t207 + t209;
t210 = t172 * t165;
t153 = atan2(-t210, -t166);
t152 = cos(t153);
t151 = sin(t153);
t199 = t151 * t210;
t138 = -t152 * t166 - t199;
t135 = 0.1e1 / t138;
t144 = 0.1e1 / t150;
t136 = 0.1e1 / t138 ^ 2;
t145 = 0.1e1 / t150 ^ 2;
t232 = t154 - 0.1e1;
t221 = t152 * t165;
t123 = (-t128 * t172 + t167) * t221 + (t233 * t166 - t195) * t151;
t231 = t123 * t135 * t136;
t183 = t166 * t209 + t207;
t212 = t167 * t174;
t196 = t165 * t212;
t132 = t183 * qJD(1) - t150 * qJD(5) + t171 * t196;
t208 = t172 * t173;
t211 = t171 * t174;
t149 = t166 * t211 - t208;
t143 = t149 ^ 2;
t142 = t143 * t145 + 0.1e1;
t224 = t145 * t149;
t189 = -qJD(1) * t166 + qJD(5);
t190 = qJD(5) * t166 - qJD(1);
t133 = -t190 * t211 + (t172 * t189 - t196) * t173;
t229 = t133 * t144 * t145;
t230 = (-t132 * t224 - t143 * t229) / t142 ^ 2;
t160 = t165 * t161;
t216 = t162 * t165;
t182 = t167 * (t160 * t162 * t163 + t216);
t217 = t161 * t172;
t187 = t205 * t217;
t228 = (t163 * t187 + t168 * t182) / t156 ^ 2;
t227 = t136 * t165;
t226 = t136 * t174;
t225 = t144 * t171;
t223 = t149 * t173;
t222 = t151 * t172;
t220 = t161 * t162;
t169 = t174 ^ 2;
t218 = t161 * t169;
t215 = t165 * t174;
t214 = t166 * t167;
t206 = qJD(1) * t172;
t131 = t136 * t218 + 0.1e1;
t204 = 0.2e1 * (-t218 * t231 + (t165 * t169 * t214 - t187) * t136) / t131 ^ 2;
t203 = 0.2e1 * t231;
t202 = -0.2e1 * t230;
t201 = t149 * t229;
t200 = t136 * t215;
t194 = 0.1e1 + t219;
t193 = t165 * t204;
t192 = -0.2e1 * t165 * t228;
t191 = t228 * t234;
t188 = t152 * t154 * t220;
t186 = t194 * t174;
t185 = t189 * t174;
t184 = t145 * t223 - t225;
t148 = -t166 * t208 + t211;
t140 = 0.1e1 / t142;
t139 = t194 * t172 * t154;
t129 = 0.1e1 / t131;
t127 = (t151 * t165 * t232 - t172 * t188) * t174;
t125 = -t166 * t222 + t221 + (t151 * t166 - t152 * t210) * t139;
t124 = -t194 * t191 + (qJD(1) * t186 + t182 * t234) * t154;
t121 = t184 * t202 * t215 + (t184 * t166 * t212 + (-t184 * t206 + ((-qJD(5) * t144 - 0.2e1 * t201) * t173 + (-t132 * t173 + (-qJD(5) * t149 + t133) * t171) * t145) * t174) * t165) * t140;
t120 = (t125 * t227 - t135 * t166) * t174 * t204 + ((-t135 * t206 + (-t125 * t167 - t123) * t226) * t166 + (-t135 * t212 - (-t124 * t152 * t172 - t233 * t151 + (t128 * t222 - t151 * t167 - t152 * t205) * t139) * t200 + (t136 * t206 + t174 * t203) * t125 - ((t124 - t205) * t151 + ((-t139 * t172 + 0.1e1) * t167 + (t139 - t172) * t128) * t152) * t166 * t226) * t165) * t129;
t1 = [t162 * t174 * t192 + (t167 * t186 - t206 * t216) * t154, t124, 0, t124, 0; (t135 * t193 + (-t135 * t214 + (qJD(1) * t127 + t123) * t227) * t129) * t172 + (t136 * t193 * t127 + (-((t192 - t214 + (t128 * t162 * t217 + t214) * t154) * t151 + (t191 * t220 - t128 * t165 + (-t160 * t198 + (t128 - 0.2e1 * t213) * t165) * t154) * t152) * t200 + (-t136 * t214 + t165 * t203) * t127 + (-t135 + ((-t168 + t169) * t188 + t232 * t199) * t136) * t165 * qJD(1)) * t129) * t174, t120, 0, t120, 0; 0.2e1 * (t144 * t183 + t148 * t224) * t230 + (0.2e1 * t148 * t201 - t190 * t144 * t208 + (t167 * t210 + t185) * t225 + (t148 * t132 + t183 * t133 - t185 * t223 - (t165 * t167 * t173 + t171 * t190) * t149 * t172) * t145) * t140, t121, 0, t121, t202 + 0.2e1 * (-t132 * t140 * t145 + (-t140 * t229 - t145 * t230) * t149) * t149;];
JaD_rot  = t1;
