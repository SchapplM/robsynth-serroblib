% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP4
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
% Datum: 2019-02-26 21:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:09:48
% EndTime: 2019-02-26 21:09:49
% DurationCPUTime: 0.73s
% Computational Cost: add. (5587->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
t172 = sin(qJ(1));
t234 = 0.2e1 * t172;
t169 = t172 ^ 2;
t167 = pkin(10) + qJ(3) + qJ(4);
t164 = sin(t167);
t160 = t164 ^ 2;
t165 = cos(t167);
t162 = 0.1e1 / t165 ^ 2;
t219 = t160 * t162;
t156 = t169 * t219 + 0.1e1;
t154 = 0.1e1 / t156;
t161 = 0.1e1 / t165;
t174 = cos(qJ(1));
t205 = qJD(1) * t174;
t195 = t164 * t205;
t168 = qJD(3) + qJD(4);
t213 = t168 * t172;
t198 = t162 * t213;
t128 = (-(-t165 * t213 - t195) * t161 + t160 * t198) * t154;
t233 = t128 - t213;
t173 = cos(qJ(5));
t207 = t173 * t174;
t171 = sin(qJ(5));
t209 = t172 * t171;
t153 = t165 * t207 + t209;
t210 = t172 * t164;
t149 = atan2(-t210, -t165);
t144 = cos(t149);
t143 = sin(t149);
t200 = t143 * t210;
t138 = -t144 * t165 - t200;
t135 = 0.1e1 / t138;
t146 = 0.1e1 / t153;
t136 = 0.1e1 / t138 ^ 2;
t147 = 0.1e1 / t153 ^ 2;
t232 = t154 - 0.1e1;
t224 = t144 * t164;
t123 = (-t128 * t172 + t168) * t224 + (t233 * t165 - t195) * t143;
t231 = t123 * t135 * t136;
t183 = t165 * t209 + t207;
t212 = t168 * t174;
t196 = t164 * t212;
t133 = t183 * qJD(1) - t153 * qJD(5) + t171 * t196;
t208 = t172 * t173;
t211 = t171 * t174;
t152 = t165 * t211 - t208;
t145 = t152 ^ 2;
t142 = t145 * t147 + 0.1e1;
t222 = t147 * t152;
t189 = -qJD(1) * t165 + qJD(5);
t190 = qJD(5) * t165 - qJD(1);
t134 = -t190 * t211 + (t172 * t189 - t196) * t173;
t228 = t134 * t146 * t147;
t230 = (-t133 * t222 - t145 * t228) / t142 ^ 2;
t159 = t164 * t160;
t216 = t161 * t164;
t182 = t168 * (t159 * t161 * t162 + t216);
t217 = t160 * t172;
t187 = t205 * t217;
t229 = (t162 * t187 + t169 * t182) / t156 ^ 2;
t227 = t136 * t164;
t226 = t136 * t174;
t225 = t143 * t172;
t223 = t146 * t171;
t221 = t152 * t173;
t220 = t160 * t161;
t170 = t174 ^ 2;
t218 = t160 * t170;
t215 = t164 * t174;
t214 = t165 * t168;
t206 = qJD(1) * t172;
t131 = t136 * t218 + 0.1e1;
t204 = 0.2e1 * (-t218 * t231 + (t164 * t170 * t214 - t187) * t136) / t131 ^ 2;
t203 = 0.2e1 * t231;
t202 = -0.2e1 * t230;
t201 = t136 * t215;
t199 = t152 * t228;
t194 = 0.1e1 + t219;
t193 = t164 * t204;
t192 = -0.2e1 * t164 * t229;
t191 = t229 * t234;
t188 = t144 * t154 * t220;
t186 = t194 * t174;
t185 = t189 * t174;
t184 = t147 * t221 - t223;
t151 = -t165 * t208 + t211;
t140 = 0.1e1 / t142;
t139 = t194 * t172 * t154;
t129 = 0.1e1 / t131;
t127 = (t143 * t164 * t232 - t172 * t188) * t174;
t125 = -t165 * t225 + t224 + (t143 * t165 - t144 * t210) * t139;
t124 = -t194 * t191 + (qJD(1) * t186 + t182 * t234) * t154;
t121 = t184 * t202 * t215 + (t184 * t165 * t212 + (-t184 * t206 + ((-qJD(5) * t146 - 0.2e1 * t199) * t173 + (-t133 * t173 + (-qJD(5) * t152 + t134) * t171) * t147) * t174) * t164) * t140;
t120 = (t125 * t227 - t135 * t165) * t174 * t204 + ((-t135 * t206 + (-t125 * t168 - t123) * t226) * t165 + (-t135 * t212 - (-t124 * t144 * t172 - t233 * t143 + (t128 * t225 - t143 * t168 - t144 * t205) * t139) * t201 + (t136 * t206 + t174 * t203) * t125 - ((t124 - t205) * t143 + ((-t139 * t172 + 0.1e1) * t168 + (t139 - t172) * t128) * t144) * t165 * t226) * t164) * t129;
t1 = [t161 * t174 * t192 + (t168 * t186 - t206 * t216) * t154, 0, t124, t124, 0, 0; (t135 * t193 + (-t135 * t214 + (qJD(1) * t127 + t123) * t227) * t129) * t172 + (t136 * t193 * t127 + (-((t192 - t214 + (t128 * t161 * t217 + t214) * t154) * t143 + (t191 * t220 - t128 * t164 + (-t159 * t198 + (t128 - 0.2e1 * t213) * t164) * t154) * t144) * t201 + (-t136 * t214 + t164 * t203) * t127 + (-t135 + ((-t169 + t170) * t188 + t232 * t200) * t136) * t164 * qJD(1)) * t129) * t174, 0, t120, t120, 0, 0; 0.2e1 * (t146 * t183 + t151 * t222) * t230 + (0.2e1 * t151 * t199 - t190 * t146 * t208 + (t168 * t210 + t185) * t223 + (t151 * t133 + t183 * t134 - t185 * t221 - (t164 * t168 * t173 + t171 * t190) * t152 * t172) * t147) * t140, 0, t121, t121, t202 + 0.2e1 * (-t133 * t147 * t140 + (-t140 * t228 - t147 * t230) * t152) * t152, 0;];
JaD_rot  = t1;
