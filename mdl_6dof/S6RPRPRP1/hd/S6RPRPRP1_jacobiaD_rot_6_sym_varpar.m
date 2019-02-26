% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:41
% EndTime: 2019-02-26 20:43:42
% DurationCPUTime: 0.75s
% Computational Cost: add. (3244->96), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->95)
t150 = qJ(3) + pkin(10);
t146 = sin(t150);
t140 = t146 ^ 2;
t148 = cos(t150);
t143 = 0.1e1 / t148 ^ 2;
t199 = t140 * t143;
t151 = qJ(1) + pkin(9);
t147 = sin(t151);
t217 = 0.2e1 * t147;
t216 = t146 * t199;
t149 = cos(t151);
t153 = cos(qJ(5));
t191 = t149 * t153;
t152 = sin(qJ(5));
t194 = t147 * t152;
t129 = t148 * t191 + t194;
t169 = qJD(5) * t148 - qJD(1);
t188 = qJD(3) * t146;
t215 = t169 * t152 + t153 * t188;
t195 = t147 * t146;
t132 = atan2(-t195, -t148);
t131 = cos(t132);
t130 = sin(t132);
t179 = t130 * t195;
t116 = -t131 * t148 - t179;
t113 = 0.1e1 / t116;
t123 = 0.1e1 / t129;
t142 = 0.1e1 / t148;
t114 = 0.1e1 / t116 ^ 2;
t124 = 0.1e1 / t129 ^ 2;
t214 = -0.2e1 * t146;
t141 = t147 ^ 2;
t135 = t141 * t199 + 0.1e1;
t133 = 0.1e1 / t135;
t213 = t133 - 0.1e1;
t189 = qJD(1) * t149;
t176 = t146 * t189;
t186 = qJD(3) * t148;
t187 = qJD(3) * t147;
t107 = (-(-t147 * t186 - t176) * t142 + t187 * t199) * t133;
t201 = t131 * t146;
t102 = (-t107 * t147 + qJD(3)) * t201 + (-t176 + (t107 - t187) * t148) * t130;
t212 = t102 * t113 * t114;
t162 = t148 * t194 + t191;
t175 = t152 * t188;
t111 = t162 * qJD(1) - t129 * qJD(5) + t149 * t175;
t192 = t149 * t152;
t193 = t147 * t153;
t128 = t148 * t192 - t193;
t122 = t128 ^ 2;
t121 = t122 * t124 + 0.1e1;
t203 = t124 * t128;
t168 = -qJD(1) * t148 + qJD(5);
t164 = t168 * t153;
t112 = t147 * t164 - t215 * t149;
t208 = t112 * t123 * t124;
t211 = (-t111 * t203 - t122 * t208) / t121 ^ 2;
t210 = t107 * t130;
t209 = t107 * t146;
t207 = t114 * t146;
t206 = t114 * t149;
t197 = t142 * t146;
t161 = qJD(3) * (t142 * t216 + t197);
t166 = t140 * t147 * t189;
t205 = (t141 * t161 + t143 * t166) / t135 ^ 2;
t173 = 0.1e1 + t199;
t118 = t173 * t147 * t133;
t204 = t118 * t147;
t202 = t130 * t148;
t200 = t140 * t142;
t145 = t149 ^ 2;
t198 = t140 * t145;
t196 = t146 * t149;
t190 = qJD(1) * t147;
t185 = qJD(3) * t149;
t110 = t114 * t198 + 0.1e1;
t184 = 0.2e1 * (-t198 * t212 + (t145 * t146 * t186 - t166) * t114) / t110 ^ 2;
t183 = 0.2e1 * t212;
t182 = -0.2e1 * t211;
t181 = t128 * t208;
t180 = t114 * t196;
t178 = t133 * t200;
t172 = t146 * t184;
t171 = t205 * t214;
t170 = t205 * t217;
t167 = t147 * t178;
t165 = t173 * t149;
t163 = -t123 * t152 + t153 * t203;
t127 = -t148 * t193 + t192;
t119 = 0.1e1 / t121;
t108 = 0.1e1 / t110;
t106 = (t213 * t146 * t130 - t131 * t167) * t149;
t104 = -t147 * t202 + t201 + (-t131 * t195 + t202) * t118;
t103 = -t173 * t170 + (qJD(1) * t165 + t161 * t217) * t133;
t1 = [t142 * t149 * t171 + (qJD(3) * t165 - t190 * t197) * t133, 0, t103, 0, 0, 0; (t113 * t172 + (-t113 * t186 + (qJD(1) * t106 + t102) * t207) * t108) * t147 + (t114 * t172 * t106 + (-((t107 * t167 + t213 * t186 + t171) * t130 + (t170 * t200 - t209 + (t209 + (t214 - t216) * t187) * t133) * t131) * t180 + (-t114 * t186 + t146 * t183) * t106 + (-t113 + ((-t141 + t145) * t131 * t178 + t213 * t179) * t114) * t146 * qJD(1)) * t108) * t149, 0 (t104 * t207 - t113 * t148) * t149 * t184 + ((-t113 * t190 + (-qJD(3) * t104 - t102) * t206) * t148 + (-t113 * t185 - (-t103 * t131 * t147 + t130 * t187 + t204 * t210 - t210 + (-qJD(3) * t130 - t131 * t189) * t118) * t180 + (t114 * t190 + t149 * t183) * t104 - ((t103 - t189) * t130 + ((0.1e1 - t204) * qJD(3) + (t118 - t147) * t107) * t131) * t148 * t206) * t146) * t108, 0, 0, 0; 0.2e1 * (t123 * t162 + t127 * t203) * t211 + (0.2e1 * t127 * t181 + (t127 * t111 + t162 * t112 + (-t215 * t147 - t149 * t164) * t128) * t124 + (t168 * t192 + (-t169 * t153 + t175) * t147) * t123) * t119, 0, t163 * t182 * t196 + (t163 * t148 * t185 + (-t163 * t190 + ((-qJD(5) * t123 - 0.2e1 * t181) * t153 + (-t111 * t153 + (-qJD(5) * t128 + t112) * t152) * t124) * t149) * t146) * t119, 0, t182 + 0.2e1 * (-t111 * t119 * t124 + (-t119 * t208 - t124 * t211) * t128) * t128, 0;];
JaD_rot  = t1;
