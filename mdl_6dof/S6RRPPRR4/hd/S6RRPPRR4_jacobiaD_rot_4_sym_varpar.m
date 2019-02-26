% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR4_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:11
% EndTime: 2019-02-26 21:30:12
% DurationCPUTime: 0.72s
% Computational Cost: add. (2643->74), mult. (7979->173), div. (436->14), fcn. (10394->11), ass. (0->83)
t160 = sin(pkin(6));
t201 = sin(pkin(11));
t202 = cos(pkin(11));
t207 = sin(qJ(2));
t208 = cos(qJ(2));
t172 = t208 * t201 + t207 * t202;
t146 = t172 * t160;
t137 = qJD(2) * t146;
t149 = t207 * t201 - t208 * t202;
t145 = t149 * t160;
t142 = 0.1e1 / t145 ^ 2;
t193 = t137 * t142;
t203 = cos(pkin(6));
t181 = t203 * t201;
t182 = t203 * t202;
t210 = t207 * t181 - t208 * t182;
t171 = t149 * qJD(2);
t161 = sin(qJ(1));
t162 = cos(qJ(1));
t129 = -t161 * t172 - t162 * t210;
t121 = atan2(t129, t145);
t116 = sin(t121);
t117 = cos(t121);
t125 = t129 ^ 2;
t120 = t125 * t142 + 0.1e1;
t118 = 0.1e1 / t120;
t141 = 0.1e1 / t145;
t195 = t129 * t141;
t209 = (t117 * t195 - t116) * t118 + t116;
t110 = t116 * t129 + t117 * t145;
t107 = 0.1e1 / t110;
t157 = 0.1e1 / t161;
t108 = 0.1e1 / t110 ^ 2;
t158 = 0.1e1 / t161 ^ 2;
t169 = t161 * t210;
t132 = -t162 * t172 + t169;
t126 = t132 ^ 2;
t106 = t126 * t108 + 0.1e1;
t147 = t208 * t181 + t207 * t182;
t140 = t147 * qJD(2);
t112 = t129 * qJD(1) - t161 * t140 - t162 * t171;
t198 = t112 * t108;
t190 = qJD(1) * t162;
t115 = qJD(1) * t169 - t162 * t140 + t161 * t171 - t172 * t190;
t176 = t115 * t141 - t129 * t193;
t101 = t176 * t118;
t180 = -t116 * t145 + t117 * t129;
t98 = t180 * t101 + t116 * t115 + t117 * t137;
t205 = t107 * t108 * t98;
t206 = 0.1e1 / t106 ^ 2 * (-t126 * t205 - t132 * t198);
t178 = -t162 * t147 + t161 * t149;
t194 = t129 * t146;
t175 = -t141 * t178 + t142 * t194;
t102 = t175 * t118;
t99 = -t180 * t102 + t116 * t178 + t117 * t146;
t204 = t132 * t99;
t192 = t141 * t193;
t200 = (t129 * t142 * t115 - t125 * t192) / t120 ^ 2;
t139 = t210 * qJD(2);
t148 = t172 * qJD(2);
t113 = t178 * qJD(1) + t161 * t139 - t162 * t148;
t177 = t161 * t147 + t162 * t149;
t127 = t177 ^ 2;
t156 = 0.1e1 / t160 ^ 2;
t124 = t127 * t158 * t156 + 0.1e1;
t159 = t157 * t158;
t199 = (-t113 * t158 * t177 - t127 * t159 * t190) * t156 / t124 ^ 2;
t197 = t116 * t132;
t196 = t117 * t132;
t191 = t158 * t162;
t189 = -0.2e1 * t206;
t188 = -0.2e1 * t205;
t187 = 0.2e1 * t200;
t186 = -0.2e1 * t141 * t200;
t179 = t162 * t139 + t161 * t148;
t155 = 0.1e1 / t160;
t138 = t160 * t171;
t122 = 0.1e1 / t124;
t114 = t177 * qJD(1) + t179;
t104 = 0.1e1 / t106;
t100 = t209 * t132;
t97 = t175 * t187 + (0.2e1 * t192 * t194 + t114 * t141 + (-t115 * t146 + t129 * t138 - t137 * t178) * t142) * t118;
t1 = [t132 * t186 + (-t112 * t141 - t132 * t193) * t118, t97, 0, 0, 0, 0; t129 * t107 * t189 + (t115 * t107 + (-t100 * t112 - t129 * t98) * t108) * t104 + ((t100 * t188 - t209 * t198) * t104 + (t100 * t189 + ((-t101 * t118 * t195 + t187) * t197 + (t129 * t186 + t101 + (-t101 + t176) * t118) * t196) * t104) * t108) * t132, 0.2e1 * (t107 * t177 - t108 * t204) * t206 + (t113 * t107 + t188 * t204 + (-t99 * t112 + t177 * t98 + (-t102 * t115 + t129 * t97 - t138 + (t102 * t145 + t178) * t101) * t196 + (t102 * t137 - t145 * t97 + t114 + (t102 * t129 - t146) * t101) * t197) * t108) * t104, 0, 0, 0, 0; (0.2e1 * (-t157 * t178 - t177 * t191) * t199 + (t179 * t157 - t113 * t191 + (-0.2e1 * t162 ^ 2 * t159 * t177 - t178 * t191) * qJD(1)) * t122) * t155 (-0.2e1 * t132 * t157 * t199 + (-t132 * t158 * t190 - t112 * t157) * t122) * t155, 0, 0, 0, 0;];
JaD_rot  = t1;
