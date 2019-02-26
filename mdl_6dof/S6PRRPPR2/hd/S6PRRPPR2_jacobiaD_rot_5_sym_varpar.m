% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPPR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:46
% EndTime: 2019-02-26 19:58:47
% DurationCPUTime: 0.79s
% Computational Cost: add. (4935->90), mult. (7372->193), div. (524->12), fcn. (9540->11), ass. (0->88)
t176 = sin(pkin(10));
t178 = cos(pkin(10));
t180 = sin(qJ(2));
t179 = cos(pkin(6));
t181 = cos(qJ(2));
t202 = t179 * t181;
t165 = -t176 * t180 + t178 * t202;
t158 = t165 * qJD(2);
t203 = t179 * t180;
t166 = t176 * t181 + t178 * t203;
t175 = qJ(3) + pkin(11);
t173 = sin(t175);
t177 = sin(pkin(6));
t206 = t177 * t178;
t194 = t173 * t206;
t174 = cos(t175);
t199 = qJD(3) * t174;
t130 = -qJD(3) * t194 + t158 * t173 + t166 * t199;
t148 = t166 * t173 + t174 * t206;
t145 = t148 ^ 2;
t205 = t177 * t180;
t156 = t173 * t205 - t174 * t179;
t154 = 0.1e1 / t156 ^ 2;
t138 = t145 * t154 + 0.1e1;
t136 = 0.1e1 / t138;
t157 = t173 * t179 + t174 * t205;
t200 = qJD(2) * t181;
t193 = t177 * t200;
t143 = qJD(3) * t157 + t173 * t193;
t153 = 0.1e1 / t156;
t210 = t148 * t154;
t118 = (-t130 * t153 + t143 * t210) * t136;
t139 = atan2(-t148, t156);
t134 = sin(t139);
t135 = cos(t139);
t191 = -t134 * t156 - t135 * t148;
t115 = t118 * t191 - t130 * t134 + t135 * t143;
t129 = -t134 * t148 + t135 * t156;
t126 = 0.1e1 / t129;
t127 = 0.1e1 / t129 ^ 2;
t219 = t115 * t126 * t127;
t195 = t176 * t203;
t168 = t178 * t181 - t195;
t207 = t176 * t177;
t189 = -t168 * t173 + t174 * t207;
t218 = -0.2e1 * t189 * t219;
t204 = t177 * t181;
t188 = -t153 * t165 + t204 * t210;
t217 = t173 * t188;
t211 = t143 * t153 * t154;
t216 = -0.2e1 * (t130 * t210 - t145 * t211) / t138 ^ 2;
t167 = t176 * t202 + t178 * t180;
t162 = 0.1e1 / t167;
t163 = 0.1e1 / t167 ^ 2;
t215 = t127 * t189;
t152 = t168 * t174 + t173 * t207;
t160 = t167 * qJD(2);
t132 = qJD(3) * t152 - t160 * t173;
t214 = t132 * t127;
t213 = t134 * t189;
t212 = t135 * t189;
t209 = t152 * t168;
t208 = t167 * t173;
t201 = qJD(2) * t180;
t146 = t189 ^ 2;
t124 = t127 * t146 + 0.1e1;
t198 = 0.2e1 * (-t146 * t219 - t189 * t214) / t124 ^ 2;
t133 = qJD(3) * t189 - t160 * t174;
t147 = t152 ^ 2;
t142 = t147 * t163 + 0.1e1;
t161 = -qJD(2) * t195 + t178 * t200;
t164 = t162 * t163;
t197 = 0.2e1 * (t133 * t152 * t163 - t147 * t161 * t164) / t142 ^ 2;
t192 = -0.2e1 * t148 * t211;
t150 = t166 * t174 - t194;
t190 = -t150 * t153 + t157 * t210;
t159 = t166 * qJD(2);
t144 = -qJD(3) * t156 + t174 * t193;
t140 = 0.1e1 / t142;
t131 = -qJD(3) * t148 + t158 * t174;
t121 = 0.1e1 / t124;
t120 = t136 * t217;
t119 = t190 * t136;
t117 = (-t134 * t165 + t135 * t204) * t173 + t191 * t120;
t116 = t119 * t191 - t134 * t150 + t135 * t157;
t114 = t190 * t216 + (t157 * t192 - t131 * t153 + (t130 * t157 + t143 * t150 + t144 * t148) * t154) * t136;
t112 = t216 * t217 + (t188 * t199 + (t192 * t204 + t153 * t159 + (t143 * t165 + (t130 * t181 - t148 * t201) * t177) * t154) * t173) * t136;
t1 = [0, t112, t114, 0, 0, 0; 0 (-t117 * t215 + t126 * t208) * t198 + ((-t161 * t173 - t167 * t199) * t126 + (-t214 + t218) * t117 + (t208 * t115 + (-t112 * t148 - t120 * t130 + (-t173 * t201 + t181 * t199) * t177 + (-t120 * t156 - t165 * t173) * t118) * t212 + (-t165 * t199 - t112 * t156 - t120 * t143 + t159 * t173 + (t120 * t148 - t173 * t204) * t118) * t213) * t127) * t121 (-t116 * t215 - t126 * t152) * t198 + (t116 * t218 + t133 * t126 + (-t152 * t115 - t116 * t132 + (-t114 * t148 - t119 * t130 + t144 + (-t119 * t156 - t150) * t118) * t212 + (-t114 * t156 - t119 * t143 - t131 + (t119 * t148 - t157) * t118) * t213) * t127) * t121, 0, 0, 0; 0 (t162 * t167 * t174 + t163 * t209) * t197 + (qJD(3) * t162 * t208 + (-t133 * t168 + t152 * t160) * t163 + (0.2e1 * t164 * t209 + (t163 * t167 - t162) * t174) * t161) * t140, -t189 * t162 * t197 + (-t161 * t163 * t189 - t132 * t162) * t140, 0, 0, 0;];
JaD_rot  = t1;
