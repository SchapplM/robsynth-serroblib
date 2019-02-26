% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:27:24
% EndTime: 2019-02-26 20:27:25
% DurationCPUTime: 0.76s
% Computational Cost: add. (3057->96), mult. (4580->216), div. (480->12), fcn. (5898->11), ass. (0->93)
t229 = sin(pkin(9));
t230 = cos(pkin(9));
t231 = sin(qJ(1));
t232 = cos(qJ(1));
t159 = t232 * t229 - t231 * t230;
t156 = t159 ^ 2;
t173 = qJ(4) + pkin(10);
t171 = sin(t173);
t167 = t171 ^ 2;
t172 = cos(t173);
t169 = 0.1e1 / t172 ^ 2;
t210 = t167 * t169;
t152 = t156 * t210 + 0.1e1;
t158 = -t229 * t231 - t230 * t232;
t154 = t158 * qJD(1);
t166 = t171 * t167;
t168 = 0.1e1 / t172;
t209 = t168 * t171;
t186 = qJD(4) * (t166 * t168 * t169 + t209);
t222 = (t154 * t159 * t210 + t156 * t186) / t152 ^ 2;
t234 = -0.2e1 * t222;
t195 = 0.1e1 + t210;
t233 = t159 * t195;
t212 = t159 * t171;
t149 = atan2(t212, t172);
t147 = sin(t149);
t148 = cos(t149);
t135 = t147 * t212 + t148 * t172;
t132 = 0.1e1 / t135;
t174 = sin(qJ(6));
t175 = cos(qJ(6));
t207 = t172 * t175;
t146 = -t158 * t207 + t159 * t174;
t140 = 0.1e1 / t146;
t133 = 0.1e1 / t135 ^ 2;
t141 = 0.1e1 / t146 ^ 2;
t157 = t158 ^ 2;
t213 = t157 * t167;
t129 = t133 * t213 + 0.1e1;
t155 = t159 * qJD(1);
t203 = qJD(4) * t172;
t150 = 0.1e1 / t152;
t197 = t159 * t203;
t205 = qJD(4) * t159;
t198 = t169 * t205;
t214 = t154 * t171;
t124 = ((t197 + t214) * t168 + t167 * t198) * t150;
t215 = t148 * t171;
t120 = (t124 * t159 - qJD(4)) * t215 + (t214 + (-t124 + t205) * t172) * t147;
t227 = t120 * t132 * t133;
t228 = (-t213 * t227 + (-t155 * t158 * t167 + t157 * t171 * t203) * t133) / t129 ^ 2;
t204 = qJD(4) * t171;
t183 = qJD(6) * t159 + t155 * t172 + t158 * t204;
t202 = qJD(6) * t172;
t189 = t158 * t202 + t154;
t125 = t174 * t183 - t175 * t189;
t208 = t172 * t174;
t145 = -t158 * t208 - t159 * t175;
t139 = t145 ^ 2;
t138 = t139 * t141 + 0.1e1;
t219 = t141 * t145;
t126 = t174 * t189 + t175 * t183;
t223 = t126 * t140 * t141;
t226 = (t125 * t219 - t139 * t223) / t138 ^ 2;
t131 = t150 * t233;
t216 = t147 * t172;
t122 = t159 * t216 - t215 + (t148 * t212 - t216) * t131;
t225 = t122 * t133;
t211 = t167 * t168;
t199 = t159 * t211;
t192 = t148 * t199;
t217 = t147 * t171;
t123 = (t217 + (t192 - t217) * t150) * t158;
t224 = t123 * t133;
t221 = t133 * t158;
t220 = t140 * t174;
t218 = t145 * t175;
t206 = t131 - t159;
t201 = -0.2e1 * t227;
t200 = 0.2e1 * t226;
t196 = t131 * t159 - 0.1e1;
t194 = -0.2e1 * t171 * t228;
t193 = 0.2e1 * t145 * t223;
t188 = t159 * t202 + t155;
t187 = t218 * t141 - t220;
t185 = t187 * t171;
t184 = qJD(6) * t158 + t154 * t172 - t159 * t204;
t144 = t158 * t174 + t159 * t207;
t143 = -t158 * t175 + t159 * t208;
t136 = 0.1e1 / t138;
t127 = 0.1e1 / t129;
t119 = t233 * t234 + (t154 * t195 + 0.2e1 * t159 * t186) * t150;
t1 = [t158 * t209 * t234 + (qJD(4) * t158 * t195 - t155 * t209) * t150, 0, 0, t119, 0, 0; t159 * t132 * t194 + (t132 * t197 + (t132 * t154 + (-t120 * t159 - t123 * t155) * t133) * t171) * t127 + (t194 * t224 + (t203 * t224 + (t123 * t201 + ((t147 * t203 + t192 * t234) * t158 + (-t147 * t155 + (t124 * t148 + 0.2e1 * t147 * t222) * t158) * t171 + (((-t124 * t199 - t203) * t158 + t155 * t171) * t147 + (-t155 * t199 + (t166 * t198 + t154 * t211 + (-t124 + 0.2e1 * t205) * t171) * t158) * t148) * t150) * t133) * t171) * t127) * t158, 0, 0, 0.2e1 * (t132 * t172 - t171 * t225) * t158 * t228 + ((t155 * t132 + (qJD(4) * t122 + t120) * t221) * t172 + (-t155 * t225 + (qJD(4) * t132 + t122 * t201 + ((t119 * t159 + t131 * t154) * t215 + (t206 * qJD(4) - t196 * t124) * t217) * t133) * t158 + ((-t119 + t154) * t147 + (t196 * qJD(4) - t206 * t124) * t148) * t172 * t221) * t171) * t127, 0, 0; (-t140 * t143 + t144 * t219) * t200 + (t144 * t193 + t188 * t140 * t175 + t184 * t220 + (t145 * t174 * t188 - t144 * t125 - t143 * t126 - t184 * t218) * t141) * t136, 0, 0, t158 * t185 * t200 + (t155 * t185 + (-t187 * t203 + ((qJD(6) * t140 + t193) * t175 + (-t125 * t175 + (qJD(6) * t145 - t126) * t174) * t141) * t171) * t158) * t136, 0, -0.2e1 * t226 + 0.2e1 * (t125 * t136 * t141 + (-t136 * t223 - t141 * t226) * t145) * t145;];
JaD_rot  = t1;
