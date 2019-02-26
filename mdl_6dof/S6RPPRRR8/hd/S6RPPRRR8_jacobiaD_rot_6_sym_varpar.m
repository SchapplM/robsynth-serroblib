% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:38:39
% EndTime: 2019-02-26 20:38:40
% DurationCPUTime: 0.70s
% Computational Cost: add. (2698->94), mult. (2734->203), div. (498->12), fcn. (3199->9), ass. (0->97)
t158 = pkin(10) + qJ(4);
t154 = sin(t158);
t150 = 0.1e1 / t154 ^ 2;
t155 = cos(t158);
t153 = t155 ^ 2;
t206 = t150 * t153;
t164 = cos(qJ(1));
t225 = 0.2e1 * t164;
t224 = t155 * t206;
t161 = t164 ^ 2;
t147 = t161 * t206 + 0.1e1;
t145 = 0.1e1 / t147;
t149 = 0.1e1 / t154;
t163 = sin(qJ(1));
t199 = qJD(1) * t163;
t188 = t155 * t199;
t195 = qJD(4) * t164;
t119 = ((t154 * t195 + t188) * t149 + t195 * t206) * t145;
t223 = -t119 + t195;
t201 = t164 * t155;
t144 = atan2(-t201, t154);
t142 = sin(t144);
t143 = cos(t144);
t128 = -t142 * t201 + t143 * t154;
t125 = 0.1e1 / t128;
t162 = qJ(5) + qJ(6);
t157 = cos(t162);
t202 = t163 * t157;
t156 = sin(t162);
t204 = t156 * t164;
t139 = t154 * t202 + t204;
t135 = 0.1e1 / t139;
t126 = 0.1e1 / t128 ^ 2;
t136 = 0.1e1 / t139 ^ 2;
t222 = t145 - 0.1e1;
t160 = t163 ^ 2;
t205 = t153 * t160;
t124 = t126 * t205 + 0.1e1;
t198 = qJD(1) * t164;
t180 = t153 * t163 * t198;
t197 = qJD(4) * t154;
t209 = t143 * t155;
t114 = (-t119 * t164 + qJD(4)) * t209 + (t154 * t223 + t188) * t142;
t220 = t114 * t125 * t126;
t221 = (-t205 * t220 + (-t155 * t160 * t197 + t180) * t126) / t124 ^ 2;
t159 = qJD(5) + qJD(6);
t183 = qJD(1) * t154 + t159;
t196 = qJD(4) * t163;
t172 = t155 * t196 + t183 * t164;
t184 = t154 * t159 + qJD(1);
t177 = t157 * t184;
t120 = t172 * t156 + t163 * t177;
t200 = t164 * t157;
t203 = t163 * t156;
t138 = t154 * t203 - t200;
t134 = t138 ^ 2;
t133 = t134 * t136 + 0.1e1;
t212 = t136 * t138;
t178 = t156 * t184;
t121 = t172 * t157 - t163 * t178;
t217 = t121 * t135 * t136;
t219 = (t120 * t212 - t134 * t217) / t133 ^ 2;
t218 = t119 * t155;
t216 = t126 * t155;
t215 = t126 * t163;
t207 = t149 * t155;
t175 = qJD(4) * (-t149 * t224 - t207);
t214 = (-t150 * t180 + t161 * t175) / t147 ^ 2;
t213 = t135 * t156;
t211 = t138 * t157;
t210 = t142 * t164;
t208 = t149 * t153;
t194 = -0.2e1 * t220;
t193 = 0.2e1 * t219;
t192 = t155 * t221;
t191 = t155 * t215;
t190 = t155 * t214;
t189 = t145 * t208;
t187 = 0.1e1 + t206;
t186 = 0.2e1 * t138 * t217;
t185 = t214 * t225;
t182 = t164 * t189;
t181 = t222 * t155 * t142;
t179 = t187 * t163;
t176 = t136 * t211 - t213;
t174 = t176 * t163;
t173 = t155 * t195 - t183 * t163;
t141 = t154 * t200 - t203;
t140 = t154 * t204 + t202;
t131 = 0.1e1 / t133;
t130 = t187 * t164 * t145;
t122 = 0.1e1 / t124;
t118 = (-t143 * t182 - t181) * t163;
t117 = t154 * t210 + t209 + (-t142 * t154 - t143 * t201) * t130;
t115 = -t187 * t185 + (-qJD(1) * t179 + t175 * t225) * t145;
t112 = -0.2e1 * t219 + 0.2e1 * (t120 * t131 * t136 + (-t131 * t217 - t136 * t219) * t138) * t138;
t1 = [-0.2e1 * t163 * t149 * t190 + (-qJD(4) * t179 + t198 * t207) * t145, 0, 0, t115, 0, 0; (0.2e1 * t125 * t192 + (t125 * t197 + (qJD(1) * t118 + t114) * t216) * t122) * t164 + (-0.2e1 * t126 * t192 * t118 + (((t119 * t182 + t222 * t197 + 0.2e1 * t190) * t142 + (t185 * t208 + t218 + (-t218 + (0.2e1 * t155 + t224) * t195) * t145) * t143) * t191 + (-t126 * t197 + t155 * t194) * t118 + (t125 + ((t160 - t161) * t143 * t189 - t164 * t181) * t126) * t155 * qJD(1)) * t122) * t163, 0, 0, 0.2e1 * (-t117 * t216 - t125 * t154) * t163 * t221 + ((t125 * t198 + (-qJD(4) * t117 - t114) * t215) * t154 + (t125 * t196 + (-t115 * t143 * t164 + t223 * t142 + (-qJD(4) * t142 + t119 * t210 + t143 * t199) * t130) * t191 + (t126 * t198 + t163 * t194) * t117 + ((-t115 - t199) * t142 + ((t130 * t164 - 0.1e1) * qJD(4) + (-t130 + t164) * t119) * t143) * t154 * t215) * t155) * t122, 0, 0; (-t135 * t140 + t141 * t212) * t193 + (t141 * t186 + t164 * t135 * t177 + t173 * t213 + (t164 * t138 * t178 - t141 * t120 - t140 * t121 - t173 * t211) * t136) * t131, 0, 0, t155 * t174 * t193 + (t174 * t197 + (-t176 * t198 + ((t135 * t159 + t186) * t157 + (-t120 * t157 + (t138 * t159 - t121) * t156) * t136) * t163) * t155) * t131, t112, t112;];
JaD_rot  = t1;
