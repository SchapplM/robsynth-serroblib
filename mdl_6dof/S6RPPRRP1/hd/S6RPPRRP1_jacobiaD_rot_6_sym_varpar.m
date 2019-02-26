% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:30:10
% EndTime: 2019-02-26 20:30:11
% DurationCPUTime: 0.74s
% Computational Cost: add. (3244->96), mult. (2519->203), div. (480->12), fcn. (2968->9), ass. (0->95)
t147 = pkin(10) + qJ(4);
t143 = sin(t147);
t137 = t143 ^ 2;
t145 = cos(t147);
t140 = 0.1e1 / t145 ^ 2;
t196 = t137 * t140;
t148 = qJ(1) + pkin(9);
t144 = sin(t148);
t214 = 0.2e1 * t144;
t213 = t143 * t196;
t146 = cos(t148);
t150 = cos(qJ(5));
t188 = t146 * t150;
t149 = sin(qJ(5));
t191 = t144 * t149;
t126 = t145 * t188 + t191;
t166 = qJD(5) * t145 - qJD(1);
t185 = qJD(4) * t143;
t212 = t166 * t149 + t150 * t185;
t192 = t144 * t143;
t129 = atan2(-t192, -t145);
t128 = cos(t129);
t127 = sin(t129);
t176 = t127 * t192;
t113 = -t128 * t145 - t176;
t110 = 0.1e1 / t113;
t120 = 0.1e1 / t126;
t139 = 0.1e1 / t145;
t111 = 0.1e1 / t113 ^ 2;
t121 = 0.1e1 / t126 ^ 2;
t211 = -0.2e1 * t143;
t138 = t144 ^ 2;
t132 = t138 * t196 + 0.1e1;
t130 = 0.1e1 / t132;
t210 = t130 - 0.1e1;
t186 = qJD(1) * t146;
t173 = t143 * t186;
t183 = qJD(4) * t145;
t184 = qJD(4) * t144;
t104 = (-(-t144 * t183 - t173) * t139 + t184 * t196) * t130;
t198 = t128 * t143;
t99 = (-t104 * t144 + qJD(4)) * t198 + (-t173 + (t104 - t184) * t145) * t127;
t209 = t110 * t111 * t99;
t159 = t145 * t191 + t188;
t172 = t149 * t185;
t108 = t159 * qJD(1) - t126 * qJD(5) + t146 * t172;
t189 = t146 * t149;
t190 = t144 * t150;
t125 = t145 * t189 - t190;
t119 = t125 ^ 2;
t118 = t119 * t121 + 0.1e1;
t200 = t121 * t125;
t165 = -qJD(1) * t145 + qJD(5);
t161 = t165 * t150;
t109 = t144 * t161 - t212 * t146;
t205 = t109 * t120 * t121;
t208 = (-t108 * t200 - t119 * t205) / t118 ^ 2;
t207 = t104 * t127;
t206 = t104 * t143;
t204 = t111 * t143;
t203 = t111 * t146;
t194 = t139 * t143;
t158 = qJD(4) * (t139 * t213 + t194);
t163 = t137 * t144 * t186;
t202 = (t138 * t158 + t140 * t163) / t132 ^ 2;
t170 = 0.1e1 + t196;
t115 = t170 * t144 * t130;
t201 = t115 * t144;
t199 = t127 * t145;
t197 = t137 * t139;
t142 = t146 ^ 2;
t195 = t137 * t142;
t193 = t143 * t146;
t187 = qJD(1) * t144;
t182 = qJD(4) * t146;
t107 = t111 * t195 + 0.1e1;
t181 = 0.2e1 / t107 ^ 2 * (-t195 * t209 + (t142 * t143 * t183 - t163) * t111);
t180 = 0.2e1 * t209;
t179 = -0.2e1 * t208;
t178 = t125 * t205;
t177 = t111 * t193;
t175 = t130 * t197;
t169 = t143 * t181;
t168 = t202 * t211;
t167 = t202 * t214;
t164 = t144 * t175;
t162 = t170 * t146;
t160 = -t120 * t149 + t150 * t200;
t124 = -t145 * t190 + t189;
t116 = 0.1e1 / t118;
t105 = 0.1e1 / t107;
t103 = (t210 * t143 * t127 - t128 * t164) * t146;
t101 = -t144 * t199 + t198 + (-t128 * t192 + t199) * t115;
t100 = -t170 * t167 + (qJD(1) * t162 + t158 * t214) * t130;
t1 = [t139 * t146 * t168 + (qJD(4) * t162 - t187 * t194) * t130, 0, 0, t100, 0, 0; (t110 * t169 + (-t110 * t183 + (qJD(1) * t103 + t99) * t204) * t105) * t144 + (t111 * t169 * t103 + (-((t104 * t164 + t210 * t183 + t168) * t127 + (t167 * t197 - t206 + (t206 + (t211 - t213) * t184) * t130) * t128) * t177 + (-t111 * t183 + t143 * t180) * t103 + (-t110 + ((-t138 + t142) * t128 * t175 + t210 * t176) * t111) * t143 * qJD(1)) * t105) * t146, 0, 0 (t101 * t204 - t110 * t145) * t146 * t181 + ((-t110 * t187 + (-qJD(4) * t101 - t99) * t203) * t145 + (-t110 * t182 - (-t100 * t128 * t144 + t127 * t184 + t201 * t207 - t207 + (-qJD(4) * t127 - t128 * t186) * t115) * t177 + (t111 * t187 + t146 * t180) * t101 - ((t100 - t186) * t127 + ((0.1e1 - t201) * qJD(4) + (t115 - t144) * t104) * t128) * t145 * t203) * t143) * t105, 0, 0; 0.2e1 * (t120 * t159 + t124 * t200) * t208 + (0.2e1 * t124 * t178 + (t124 * t108 + t159 * t109 + (-t212 * t144 - t146 * t161) * t125) * t121 + (t165 * t189 + (-t166 * t150 + t172) * t144) * t120) * t116, 0, 0, t160 * t179 * t193 + (t160 * t145 * t182 + (-t160 * t187 + ((-qJD(5) * t120 - 0.2e1 * t178) * t150 + (-t108 * t150 + (-qJD(5) * t125 + t109) * t149) * t121) * t146) * t143) * t116, t179 + 0.2e1 * (-t108 * t116 * t121 + (-t116 * t205 - t121 * t208) * t125) * t125, 0;];
JaD_rot  = t1;
