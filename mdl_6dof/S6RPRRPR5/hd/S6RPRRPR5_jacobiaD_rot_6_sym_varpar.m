% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:13
% EndTime: 2019-02-26 21:03:14
% DurationCPUTime: 0.75s
% Computational Cost: add. (5302->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->95)
t174 = sin(qJ(1));
t169 = pkin(10) + qJ(3) + qJ(4);
t167 = sin(t169);
t163 = 0.1e1 / t167 ^ 2;
t168 = cos(t169);
t166 = t168 ^ 2;
t220 = t163 * t166;
t195 = 0.1e1 + t220;
t235 = t174 * t195;
t171 = t174 ^ 2;
t160 = t171 * t220 + 0.1e1;
t158 = 0.1e1 / t160;
t162 = 0.1e1 / t167;
t176 = cos(qJ(1));
t207 = qJD(1) * t176;
t196 = t168 * t207;
t170 = qJD(3) + qJD(4);
t215 = t170 * t174;
t198 = t163 * t215;
t132 = ((t167 * t215 - t196) * t162 + t166 * t198) * t158;
t234 = -t132 + t215;
t191 = qJD(1) * t167 + qJD(6);
t214 = t170 * t176;
t233 = -t168 * t214 + t191 * t174;
t213 = t174 * t168;
t157 = atan2(-t213, t167);
t148 = cos(t157);
t147 = sin(t157);
t200 = t147 * t213;
t142 = t148 * t167 - t200;
t139 = 0.1e1 / t142;
t173 = sin(qJ(6));
t210 = t176 * t173;
t175 = cos(qJ(6));
t211 = t174 * t175;
t154 = t167 * t210 + t211;
t150 = 0.1e1 / t154;
t140 = 0.1e1 / t142 ^ 2;
t151 = 0.1e1 / t154 ^ 2;
t232 = t158 - 0.1e1;
t224 = t148 * t168;
t127 = (-t132 * t174 + t170) * t224 + (t234 * t167 - t196) * t147;
t231 = t127 * t139 * t140;
t192 = qJD(6) * t167 + qJD(1);
t187 = t192 * t176;
t137 = t173 * t187 + t233 * t175;
t209 = t176 * t175;
t212 = t174 * t173;
t153 = -t167 * t209 + t212;
t149 = t153 ^ 2;
t146 = t149 * t151 + 0.1e1;
t223 = t151 * t153;
t138 = -t233 * t173 + t175 * t187;
t228 = t138 * t150 * t151;
t230 = (t137 * t223 - t149 * t228) / t146 ^ 2;
t165 = t168 * t166;
t221 = t162 * t168;
t185 = t170 * (-t162 * t163 * t165 - t221);
t218 = t166 * t174;
t189 = t207 * t218;
t229 = (t163 * t189 + t171 * t185) / t160 ^ 2;
t227 = t140 * t168;
t226 = t140 * t176;
t225 = t147 * t174;
t222 = t153 * t173;
t172 = t176 ^ 2;
t219 = t166 * t172;
t217 = t167 * t170;
t216 = t168 * t170;
t208 = qJD(1) * t174;
t135 = t140 * t219 + 0.1e1;
t206 = 0.2e1 * (-t219 * t231 + (-t167 * t172 * t216 - t189) * t140) / t135 ^ 2;
t205 = 0.2e1 * t231;
t204 = 0.2e1 * t230;
t203 = -0.2e1 * t229;
t202 = t168 * t229;
t201 = t168 * t226;
t199 = t162 * t218;
t194 = t168 * t206;
t193 = 0.2e1 * t153 * t228;
t190 = t148 * t158 * t162 * t166;
t188 = t195 * t176;
t186 = t150 * t175 + t151 * t222;
t184 = t186 * t176;
t156 = -t167 * t212 + t209;
t155 = t167 * t211 + t210;
t144 = 0.1e1 / t146;
t143 = t158 * t235;
t133 = 0.1e1 / t135;
t131 = (t232 * t168 * t147 + t174 * t190) * t176;
t129 = t167 * t225 + t224 + (-t147 * t167 - t148 * t213) * t143;
t128 = t203 * t235 + (qJD(1) * t188 + 0.2e1 * t174 * t185) * t158;
t125 = t168 * t184 * t204 + (t184 * t217 + (t186 * t208 + ((qJD(6) * t150 + t193) * t173 + (-t137 * t173 + (-qJD(6) * t153 + t138) * t175) * t151) * t176) * t168) * t144;
t124 = (t129 * t227 + t139 * t167) * t176 * t206 + ((t139 * t208 + (t129 * t170 + t127) * t226) * t167 + (-t139 * t214 - (-t128 * t148 * t174 + t234 * t147 + (t132 * t225 - t147 * t170 - t148 * t207) * t143) * t201 + (t140 * t208 + t176 * t205) * t129 - ((-t128 + t207) * t147 + ((t143 * t174 - 0.1e1) * t170 + (-t143 + t174) * t132) * t148) * t167 * t226) * t168) * t133;
t1 = [0.2e1 * t176 * t162 * t202 + (t170 * t188 + t208 * t221) * t158, 0, t128, t128, 0, 0; (t139 * t194 + (t139 * t217 + (qJD(1) * t131 + t127) * t227) * t133) * t174 + (t140 * t194 * t131 + (-((-0.2e1 * t202 + t217 + (-t132 * t199 - t217) * t158) * t147 + (t199 * t203 - t132 * t168 + (-t165 * t198 + (t132 - 0.2e1 * t215) * t168) * t158) * t148) * t201 + (t140 * t217 + t168 * t205) * t131 + (-t139 + ((t171 - t172) * t190 + t232 * t200) * t140) * t168 * qJD(1)) * t133) * t176, 0, t124, t124, 0, 0; (-t150 * t155 + t156 * t223) * t204 + (t156 * t193 + (-t156 * t137 - t155 * t138 + t192 * t153 * t211 - (-t170 * t213 - t191 * t176) * t222) * t151 + (t191 * t209 + (-t192 * t173 + t175 * t216) * t174) * t150) * t144, 0, t125, t125, 0, -0.2e1 * t230 + 0.2e1 * (t137 * t151 * t144 + (-t144 * t228 - t151 * t230) * t153) * t153;];
JaD_rot  = t1;
