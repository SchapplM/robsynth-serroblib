% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP5
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
% Datum: 2019-02-26 21:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:10:24
% EndTime: 2019-02-26 21:10:24
% DurationCPUTime: 0.75s
% Computational Cost: add. (5587->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
t174 = sin(qJ(1));
t236 = 0.2e1 * t174;
t171 = t174 ^ 2;
t169 = pkin(10) + qJ(3) + qJ(4);
t166 = sin(t169);
t162 = t166 ^ 2;
t167 = cos(t169);
t164 = 0.1e1 / t167 ^ 2;
t221 = t162 * t164;
t158 = t171 * t221 + 0.1e1;
t156 = 0.1e1 / t158;
t163 = 0.1e1 / t167;
t176 = cos(qJ(1));
t207 = qJD(1) * t176;
t197 = t166 * t207;
t170 = qJD(3) + qJD(4);
t215 = t170 * t174;
t200 = t164 * t215;
t130 = (-(-t167 * t215 - t197) * t163 + t162 * t200) * t156;
t235 = t130 - t215;
t175 = cos(qJ(5));
t209 = t175 * t176;
t173 = sin(qJ(5));
t211 = t174 * t173;
t155 = t167 * t209 + t211;
t212 = t174 * t166;
t151 = atan2(-t212, -t167);
t146 = cos(t151);
t145 = sin(t151);
t202 = t145 * t212;
t140 = -t146 * t167 - t202;
t137 = 0.1e1 / t140;
t148 = 0.1e1 / t155;
t138 = 0.1e1 / t140 ^ 2;
t149 = 0.1e1 / t155 ^ 2;
t234 = t156 - 0.1e1;
t226 = t146 * t166;
t125 = (-t130 * t174 + t170) * t226 + (t235 * t167 - t197) * t145;
t233 = t125 * t137 * t138;
t185 = t167 * t211 + t209;
t214 = t170 * t176;
t198 = t166 * t214;
t135 = t185 * qJD(1) - t155 * qJD(5) + t173 * t198;
t210 = t174 * t175;
t213 = t173 * t176;
t154 = t167 * t213 - t210;
t147 = t154 ^ 2;
t144 = t147 * t149 + 0.1e1;
t224 = t149 * t154;
t191 = -qJD(1) * t167 + qJD(5);
t192 = qJD(5) * t167 - qJD(1);
t136 = -t192 * t213 + (t174 * t191 - t198) * t175;
t230 = t136 * t148 * t149;
t232 = (-t135 * t224 - t147 * t230) / t144 ^ 2;
t161 = t166 * t162;
t218 = t163 * t166;
t184 = t170 * (t161 * t163 * t164 + t218);
t219 = t162 * t174;
t189 = t207 * t219;
t231 = (t164 * t189 + t171 * t184) / t158 ^ 2;
t229 = t138 * t166;
t228 = t138 * t176;
t227 = t145 * t174;
t225 = t148 * t173;
t223 = t154 * t175;
t222 = t162 * t163;
t172 = t176 ^ 2;
t220 = t162 * t172;
t217 = t166 * t176;
t216 = t167 * t170;
t208 = qJD(1) * t174;
t133 = t138 * t220 + 0.1e1;
t206 = 0.2e1 * (-t220 * t233 + (t166 * t172 * t216 - t189) * t138) / t133 ^ 2;
t205 = 0.2e1 * t233;
t204 = -0.2e1 * t232;
t203 = t138 * t217;
t201 = t154 * t230;
t196 = 0.1e1 + t221;
t195 = t166 * t206;
t194 = -0.2e1 * t166 * t231;
t193 = t231 * t236;
t190 = t146 * t156 * t222;
t188 = t196 * t176;
t187 = t191 * t176;
t186 = t149 * t223 - t225;
t153 = -t167 * t210 + t213;
t142 = 0.1e1 / t144;
t141 = t196 * t174 * t156;
t131 = 0.1e1 / t133;
t129 = (t145 * t166 * t234 - t174 * t190) * t176;
t127 = -t167 * t227 + t226 + (t145 * t167 - t146 * t212) * t141;
t126 = -t196 * t193 + (qJD(1) * t188 + t184 * t236) * t156;
t123 = t186 * t204 * t217 + (t186 * t167 * t214 + (-t186 * t208 + ((-qJD(5) * t148 - 0.2e1 * t201) * t175 + (-t135 * t175 + (-qJD(5) * t154 + t136) * t173) * t149) * t176) * t166) * t142;
t122 = (t127 * t229 - t137 * t167) * t176 * t206 + ((-t137 * t208 + (-t127 * t170 - t125) * t228) * t167 + (-t137 * t214 - (-t126 * t146 * t174 - t235 * t145 + (t130 * t227 - t145 * t170 - t146 * t207) * t141) * t203 + (t138 * t208 + t176 * t205) * t127 - ((t126 - t207) * t145 + ((-t141 * t174 + 0.1e1) * t170 + (t141 - t174) * t130) * t146) * t167 * t228) * t166) * t131;
t1 = [t163 * t176 * t194 + (t170 * t188 - t208 * t218) * t156, 0, t126, t126, 0, 0; (t137 * t195 + (-t137 * t216 + (qJD(1) * t129 + t125) * t229) * t131) * t174 + (t138 * t195 * t129 + (-((t194 - t216 + (t130 * t163 * t219 + t216) * t156) * t145 + (t193 * t222 - t130 * t166 + (-t161 * t200 + (t130 - 0.2e1 * t215) * t166) * t156) * t146) * t203 + (-t138 * t216 + t166 * t205) * t129 + (-t137 + ((-t171 + t172) * t190 + t234 * t202) * t138) * t166 * qJD(1)) * t131) * t176, 0, t122, t122, 0, 0; 0.2e1 * (t148 * t185 + t153 * t224) * t232 + (0.2e1 * t153 * t201 - t192 * t148 * t210 + (t170 * t212 + t187) * t225 + (t153 * t135 + t185 * t136 - t187 * t223 - (t166 * t170 * t175 + t173 * t192) * t154 * t174) * t149) * t142, 0, t123, t123, t204 + 0.2e1 * (-t135 * t149 * t142 + (-t142 * t230 - t149 * t232) * t154) * t154, 0;];
JaD_rot  = t1;
