% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRP3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:36
% EndTime: 2019-02-26 22:10:36
% DurationCPUTime: 0.82s
% Computational Cost: add. (4131->97), mult. (3810->205), div. (753->12), fcn. (4455->9), ass. (0->96)
t170 = sin(qJ(1));
t231 = 0.2e1 * t170;
t167 = t170 ^ 2;
t169 = qJ(2) + qJ(3);
t163 = sin(t169);
t157 = t163 ^ 2;
t164 = cos(t169);
t159 = 0.1e1 / t164 ^ 2;
t216 = t157 * t159;
t153 = t167 * t216 + 0.1e1;
t150 = 0.1e1 / t153;
t158 = 0.1e1 / t164;
t171 = cos(qJ(1));
t202 = qJD(1) * t171;
t192 = t163 * t202;
t166 = qJD(2) + qJD(3);
t210 = t166 * t170;
t195 = t159 * t210;
t124 = (-(-t164 * t210 - t192) * t158 + t157 * t195) * t150;
t230 = t124 - t210;
t165 = pkin(10) + qJ(5);
t162 = cos(t165);
t204 = t171 * t162;
t161 = sin(t165);
t208 = t170 * t161;
t146 = t164 * t204 + t208;
t206 = t170 * t163;
t149 = atan2(-t206, -t164);
t148 = cos(t149);
t147 = sin(t149);
t196 = t147 * t206;
t134 = -t148 * t164 - t196;
t131 = 0.1e1 / t134;
t140 = 0.1e1 / t146;
t132 = 0.1e1 / t134 ^ 2;
t141 = 0.1e1 / t146 ^ 2;
t229 = t150 - 0.1e1;
t218 = t148 * t163;
t119 = (-t124 * t170 + t166) * t218 + (t230 * t164 - t192) * t147;
t228 = t119 * t131 * t132;
t181 = t164 * t208 + t204;
t209 = t166 * t171;
t193 = t163 * t209;
t125 = t181 * qJD(1) - t146 * qJD(5) + t161 * t193;
t205 = t171 * t161;
t207 = t170 * t162;
t145 = t164 * t205 - t207;
t139 = t145 ^ 2;
t137 = t139 * t141 + 0.1e1;
t221 = t141 * t145;
t186 = -qJD(1) * t164 + qJD(5);
t187 = qJD(5) * t164 - qJD(1);
t126 = -t187 * t205 + (t186 * t170 - t193) * t162;
t226 = t126 * t140 * t141;
t227 = (-t125 * t221 - t139 * t226) / t137 ^ 2;
t156 = t163 * t157;
t213 = t158 * t163;
t180 = t166 * (t156 * t158 * t159 + t213);
t214 = t157 * t170;
t184 = t202 * t214;
t225 = (t159 * t184 + t167 * t180) / t153 ^ 2;
t224 = t132 * t163;
t223 = t132 * t171;
t222 = t140 * t161;
t220 = t145 * t162;
t219 = t147 * t170;
t217 = t157 * t158;
t168 = t171 ^ 2;
t215 = t157 * t168;
t212 = t163 * t171;
t211 = t164 * t166;
t203 = qJD(1) * t170;
t129 = t132 * t215 + 0.1e1;
t201 = 0.2e1 * (-t215 * t228 + (t163 * t168 * t211 - t184) * t132) / t129 ^ 2;
t200 = 0.2e1 * t228;
t199 = -0.2e1 * t227;
t198 = t132 * t212;
t197 = t145 * t226;
t191 = 0.1e1 + t216;
t190 = t163 * t201;
t189 = -0.2e1 * t163 * t225;
t188 = t225 * t231;
t185 = t148 * t150 * t217;
t183 = t191 * t171;
t182 = t141 * t220 - t222;
t179 = t166 * t206 + t186 * t171;
t144 = -t164 * t207 + t205;
t138 = t191 * t170 * t150;
t135 = 0.1e1 / t137;
t127 = 0.1e1 / t129;
t123 = (t229 * t163 * t147 - t170 * t185) * t171;
t122 = -t164 * t219 + t218 + (t147 * t164 - t148 * t206) * t138;
t120 = -t191 * t188 + (qJD(1) * t183 + t180 * t231) * t150;
t117 = t182 * t199 * t212 + (t182 * t164 * t209 + (-t182 * t203 + ((-qJD(5) * t140 - 0.2e1 * t197) * t162 + (-t125 * t162 + (-qJD(5) * t145 + t126) * t161) * t141) * t171) * t163) * t135;
t116 = (t122 * t224 - t131 * t164) * t171 * t201 + ((-t131 * t203 + (-t122 * t166 - t119) * t223) * t164 + (-t131 * t209 - (-t120 * t148 * t170 - t230 * t147 + (t124 * t219 - t147 * t166 - t148 * t202) * t138) * t198 + (t132 * t203 + t171 * t200) * t122 - ((t120 - t202) * t147 + ((-t138 * t170 + 0.1e1) * t166 + (t138 - t170) * t124) * t148) * t164 * t223) * t163) * t127;
t1 = [t171 * t158 * t189 + (t166 * t183 - t203 * t213) * t150, t120, t120, 0, 0, 0; (t131 * t190 + (-t131 * t211 + (qJD(1) * t123 + t119) * t224) * t127) * t170 + (t132 * t190 * t123 + (-((t189 - t211 + (t124 * t158 * t214 + t211) * t150) * t147 + (t188 * t217 - t124 * t163 + (-t156 * t195 + (t124 - 0.2e1 * t210) * t163) * t150) * t148) * t198 + (-t132 * t211 + t163 * t200) * t123 + (-t131 + ((-t167 + t168) * t185 + t229 * t196) * t132) * t163 * qJD(1)) * t127) * t171, t116, t116, 0, 0, 0; 0.2e1 * (t140 * t181 + t144 * t221) * t227 + (0.2e1 * t144 * t197 - t187 * t140 * t207 + t179 * t222 + (-t187 * t145 * t208 + t144 * t125 + t126 * t181 - t179 * t220) * t141) * t135, t117, t117, 0, t199 + 0.2e1 * (-t125 * t141 * t135 + (-t135 * t226 - t141 * t227) * t145) * t145, 0;];
JaD_rot  = t1;
