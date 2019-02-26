% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:48
% EndTime: 2019-02-26 22:41:49
% DurationCPUTime: 0.75s
% Computational Cost: add. (2348->96), mult. (2949->204), div. (516->12), fcn. (3430->9), ass. (0->95)
t161 = sin(qJ(2));
t155 = t161 ^ 2;
t163 = cos(qJ(2));
t158 = 0.1e1 / t163 ^ 2;
t208 = t155 * t158;
t162 = sin(qJ(1));
t226 = 0.2e1 * t162;
t225 = t161 * t208;
t153 = qJ(3) + qJ(4) + qJ(5);
t151 = cos(t153);
t164 = cos(qJ(1));
t200 = t163 * t164;
t150 = sin(t153);
t204 = t162 * t150;
t140 = t151 * t200 + t204;
t202 = t162 * t161;
t145 = atan2(-t202, -t163);
t144 = cos(t145);
t143 = sin(t145);
t189 = t143 * t202;
t130 = -t144 * t163 - t189;
t127 = 0.1e1 / t130;
t134 = 0.1e1 / t140;
t157 = 0.1e1 / t163;
t128 = 0.1e1 / t130 ^ 2;
t135 = 0.1e1 / t140 ^ 2;
t224 = -0.2e1 * t161;
t156 = t162 ^ 2;
t148 = t156 * t208 + 0.1e1;
t146 = 0.1e1 / t148;
t223 = t146 - 0.1e1;
t152 = qJD(3) + qJD(4) + qJD(5);
t201 = t162 * t163;
t174 = t150 * t201 + t151 * t164;
t195 = qJD(2) * t164;
t185 = t161 * t195;
t118 = t174 * qJD(1) - t140 * t152 + t150 * t185;
t203 = t162 * t151;
t139 = t150 * t200 - t203;
t133 = t139 ^ 2;
t123 = t133 * t135 + 0.1e1;
t213 = t135 * t139;
t179 = -qJD(1) * t163 + t152;
t180 = t152 * t163 - qJD(1);
t210 = t150 * t164;
t119 = -t180 * t210 + (t179 * t162 - t185) * t151;
t220 = t119 * t134 * t135;
t222 = (-t118 * t213 - t133 * t220) / t123 ^ 2;
t198 = qJD(1) * t164;
t186 = t161 * t198;
t196 = qJD(2) * t163;
t197 = qJD(2) * t162;
t120 = (-(-t162 * t196 - t186) * t157 + t197 * t208) * t146;
t211 = t144 * t161;
t114 = (-t120 * t162 + qJD(2)) * t211 + (-t186 + (t120 - t197) * t163) * t143;
t221 = t114 * t127 * t128;
t219 = t120 * t143;
t218 = t120 * t161;
t217 = t128 * t161;
t206 = t157 * t161;
t173 = qJD(2) * (t157 * t225 + t206);
t177 = t155 * t162 * t198;
t216 = (t156 * t173 + t158 * t177) / t148 ^ 2;
t184 = 0.1e1 + t208;
t132 = t184 * t162 * t146;
t215 = t132 * t162;
t214 = t134 * t150;
t212 = t139 * t151;
t209 = t155 * t157;
t160 = t164 ^ 2;
t207 = t155 * t160;
t205 = t161 * t164;
t199 = qJD(1) * t162;
t126 = t128 * t207 + 0.1e1;
t194 = 0.2e1 * (-t207 * t221 + (t160 * t161 * t196 - t177) * t128) / t126 ^ 2;
t193 = -0.2e1 * t222;
t192 = 0.2e1 * t221;
t191 = t139 * t220;
t190 = t128 * t205;
t188 = t146 * t209;
t183 = t161 * t194;
t182 = t216 * t224;
t181 = t216 * t226;
t178 = t162 * t188;
t176 = t184 * t164;
t175 = t135 * t212 - t214;
t172 = t161 * t197 + t179 * t164;
t138 = -t151 * t201 + t210;
t124 = 0.1e1 / t126;
t121 = 0.1e1 / t123;
t117 = (t223 * t161 * t143 - t144 * t178) * t164;
t116 = -t143 * t201 + t211 + (t143 * t163 - t144 * t202) * t132;
t115 = -t184 * t181 + (qJD(1) * t176 + t173 * t226) * t146;
t111 = t193 + 0.2e1 * (-t118 * t121 * t135 + (-t121 * t220 - t135 * t222) * t139) * t139;
t1 = [t157 * t164 * t182 + (qJD(2) * t176 - t199 * t206) * t146, t115, 0, 0, 0, 0; (t127 * t183 + (-t127 * t196 + (qJD(1) * t117 + t114) * t217) * t124) * t162 + (t128 * t183 * t117 + (-((t120 * t178 + t223 * t196 + t182) * t143 + (t181 * t209 - t218 + (t218 + (t224 - t225) * t197) * t146) * t144) * t190 + (-t128 * t196 + t161 * t192) * t117 + (-t127 + ((-t156 + t160) * t144 * t188 + t223 * t189) * t128) * t161 * qJD(1)) * t124) * t164 (t116 * t217 - t127 * t163) * t164 * t194 + ((-t127 * t199 + (-qJD(2) * t116 - t114) * t164 * t128) * t163 + (-t127 * t195 - (-t115 * t144 * t162 + t143 * t197 + t215 * t219 - t219 + (-qJD(2) * t143 - t144 * t198) * t132) * t190 + (t128 * t199 + t164 * t192) * t116 - ((t115 - t198) * t143 + ((0.1e1 - t215) * qJD(2) + (t132 - t162) * t120) * t144) * t128 * t200) * t161) * t124, 0, 0, 0, 0; 0.2e1 * (t134 * t174 + t138 * t213) * t222 + (0.2e1 * t138 * t191 - t180 * t134 * t203 + t172 * t214 + (-t180 * t139 * t204 + t138 * t118 + t119 * t174 - t172 * t212) * t135) * t121, t175 * t193 * t205 + (t175 * t163 * t195 + (-t175 * t199 + ((-t134 * t152 - 0.2e1 * t191) * t151 + (-t118 * t151 + (-t139 * t152 + t119) * t150) * t135) * t164) * t161) * t121, t111, t111, t111, 0;];
JaD_rot  = t1;
