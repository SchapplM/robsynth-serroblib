% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:10
% EndTime: 2019-02-26 21:47:10
% DurationCPUTime: 0.74s
% Computational Cost: add. (2887->96), mult. (2734->202), div. (498->12), fcn. (3199->9), ass. (0->96)
t166 = qJ(2) + pkin(10);
t161 = sin(t166);
t157 = t161 ^ 2;
t162 = cos(t166);
t159 = 0.1e1 / t162 ^ 2;
t215 = t157 * t159;
t170 = sin(qJ(1));
t233 = 0.2e1 * t170;
t232 = t161 * t215;
t167 = t170 ^ 2;
t153 = t167 * t215 + 0.1e1;
t151 = 0.1e1 / t153;
t158 = 0.1e1 / t162;
t171 = cos(qJ(1));
t205 = qJD(1) * t171;
t193 = t161 * t205;
t203 = qJD(2) * t170;
t124 = (-(-t162 * t203 - t193) * t158 + t203 * t215) * t151;
t231 = t124 - t203;
t169 = qJ(4) + qJ(5);
t163 = sin(t169);
t208 = t170 * t163;
t164 = cos(t169);
t210 = t164 * t171;
t146 = t162 * t210 + t208;
t209 = t170 * t161;
t149 = atan2(-t209, -t162);
t148 = cos(t149);
t147 = sin(t149);
t196 = t147 * t209;
t133 = -t148 * t162 - t196;
t130 = 0.1e1 / t133;
t140 = 0.1e1 / t146;
t131 = 0.1e1 / t133 ^ 2;
t141 = 0.1e1 / t146 ^ 2;
t230 = -0.2e1 * t161;
t229 = t151 - 0.1e1;
t217 = t148 * t161;
t119 = (-t124 * t170 + qJD(2)) * t217 + (t162 * t231 - t193) * t147;
t228 = t119 * t130 * t131;
t165 = qJD(4) + qJD(5);
t181 = t162 * t208 + t210;
t202 = qJD(2) * t171;
t192 = t161 * t202;
t125 = t181 * qJD(1) - t146 * t165 + t163 * t192;
t207 = t170 * t164;
t211 = t163 * t171;
t145 = t162 * t211 - t207;
t139 = t145 ^ 2;
t138 = t139 * t141 + 0.1e1;
t220 = t141 * t145;
t186 = -qJD(1) * t162 + t165;
t187 = t162 * t165 - qJD(1);
t126 = -t187 * t211 + (t186 * t170 - t192) * t164;
t225 = t126 * t140 * t141;
t227 = (-t125 * t220 - t139 * t225) / t138 ^ 2;
t226 = t124 * t161;
t224 = t131 * t161;
t223 = t131 * t171;
t213 = t158 * t161;
t180 = qJD(2) * (t158 * t232 + t213);
t184 = t157 * t170 * t205;
t222 = (t159 * t184 + t167 * t180) / t153 ^ 2;
t221 = t140 * t163;
t219 = t145 * t164;
t218 = t147 * t170;
t216 = t157 * t158;
t168 = t171 ^ 2;
t214 = t157 * t168;
t212 = t161 * t171;
t206 = qJD(1) * t170;
t204 = qJD(2) * t162;
t129 = t131 * t214 + 0.1e1;
t201 = 0.2e1 * (-t214 * t228 + (t161 * t168 * t204 - t184) * t131) / t129 ^ 2;
t200 = 0.2e1 * t228;
t199 = -0.2e1 * t227;
t198 = t145 * t225;
t197 = t131 * t212;
t195 = t151 * t216;
t191 = 0.1e1 + t215;
t190 = t161 * t201;
t189 = t222 * t230;
t188 = t222 * t233;
t185 = t170 * t195;
t183 = t191 * t171;
t182 = t141 * t219 - t221;
t179 = t161 * t203 + t186 * t171;
t144 = -t162 * t207 + t211;
t136 = 0.1e1 / t138;
t135 = t191 * t170 * t151;
t127 = 0.1e1 / t129;
t123 = (t229 * t161 * t147 - t148 * t185) * t171;
t122 = -t162 * t218 + t217 + (t147 * t162 - t148 * t209) * t135;
t120 = -t191 * t188 + (qJD(1) * t183 + t180 * t233) * t151;
t117 = t199 + 0.2e1 * (-t125 * t136 * t141 + (-t136 * t225 - t141 * t227) * t145) * t145;
t1 = [t158 * t171 * t189 + (qJD(2) * t183 - t206 * t213) * t151, t120, 0, 0, 0, 0; (t130 * t190 + (-t130 * t204 + (qJD(1) * t123 + t119) * t224) * t127) * t170 + (t131 * t190 * t123 + (-((t124 * t185 + t229 * t204 + t189) * t147 + (t188 * t216 - t226 + (t226 + (t230 - t232) * t203) * t151) * t148) * t197 + (-t131 * t204 + t161 * t200) * t123 + (-t130 + ((-t167 + t168) * t148 * t195 + t229 * t196) * t131) * t161 * qJD(1)) * t127) * t171 (t122 * t224 - t130 * t162) * t171 * t201 + ((-t130 * t206 + (-qJD(2) * t122 - t119) * t223) * t162 + (-t130 * t202 - (-t120 * t148 * t170 - t231 * t147 + (-qJD(2) * t147 + t124 * t218 - t148 * t205) * t135) * t197 + (t131 * t206 + t171 * t200) * t122 - ((t120 - t205) * t147 + ((-t135 * t170 + 0.1e1) * qJD(2) + (t135 - t170) * t124) * t148) * t162 * t223) * t161) * t127, 0, 0, 0, 0; 0.2e1 * (t140 * t181 + t144 * t220) * t227 + (0.2e1 * t144 * t198 - t187 * t140 * t207 + t179 * t221 + (-t187 * t145 * t208 + t144 * t125 + t126 * t181 - t179 * t219) * t141) * t136, t182 * t199 * t212 + (t182 * t162 * t202 + (-t182 * t206 + ((-t140 * t165 - 0.2e1 * t198) * t164 + (-t125 * t164 + (-t145 * t165 + t126) * t163) * t141) * t171) * t161) * t136, 0, t117, t117, 0;];
JaD_rot  = t1;
