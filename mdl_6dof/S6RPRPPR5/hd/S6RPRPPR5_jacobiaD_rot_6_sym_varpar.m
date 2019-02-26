% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:29
% EndTime: 2019-02-26 20:41:30
% DurationCPUTime: 0.67s
% Computational Cost: add. (2429->94), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->95)
t143 = pkin(9) + qJ(3);
t139 = sin(t143);
t134 = 0.1e1 / t139 ^ 2;
t141 = cos(t143);
t137 = t141 ^ 2;
t191 = t134 * t137;
t210 = t141 * t191;
t146 = sin(qJ(1));
t167 = 0.1e1 + t191;
t209 = t146 * t167;
t163 = qJD(1) * t139 + qJD(6);
t147 = cos(qJ(1));
t179 = qJD(3) * t147;
t208 = -t141 * t179 + t163 * t146;
t180 = qJD(3) * t146;
t207 = t141 * t180 + t163 * t147;
t184 = t146 * t141;
t128 = atan2(-t184, t139);
t127 = cos(t128);
t126 = sin(t128);
t172 = t126 * t184;
t112 = t127 * t139 - t172;
t109 = 0.1e1 / t112;
t142 = pkin(10) + qJ(6);
t138 = sin(t142);
t140 = cos(t142);
t185 = t146 * t140;
t187 = t139 * t147;
t123 = t138 * t187 + t185;
t119 = 0.1e1 / t123;
t133 = 0.1e1 / t139;
t110 = 0.1e1 / t112 ^ 2;
t120 = 0.1e1 / t123 ^ 2;
t144 = t146 ^ 2;
t131 = t144 * t191 + 0.1e1;
t129 = 0.1e1 / t131;
t206 = t129 - 0.1e1;
t182 = qJD(1) * t147;
t170 = t141 * t182;
t103 = ((t139 * t180 - t170) * t133 + t180 * t191) * t129;
t193 = t127 * t141;
t98 = (-t103 * t146 + qJD(3)) * t193 + (-t170 + (-t103 + t180) * t139) * t126;
t205 = t109 * t110 * t98;
t164 = qJD(6) * t139 + qJD(1);
t158 = t164 * t147;
t104 = t138 * t158 + t140 * t208;
t186 = t140 * t147;
t122 = t138 * t146 - t139 * t186;
t118 = t122 ^ 2;
t117 = t118 * t120 + 0.1e1;
t195 = t120 * t122;
t105 = -t138 * t208 + t140 * t158;
t201 = t105 * t119 * t120;
t204 = (t104 * t195 - t118 * t201) / t117 ^ 2;
t203 = t103 * t126;
t202 = t103 * t141;
t200 = t110 * t141;
t199 = t110 * t147;
t192 = t133 * t141;
t156 = qJD(3) * (-t133 * t210 - t192);
t189 = t137 * t146;
t161 = t182 * t189;
t198 = (t134 * t161 + t144 * t156) / t131 ^ 2;
t116 = t129 * t209;
t197 = t116 * t146;
t196 = t119 * t140;
t194 = t122 * t138;
t145 = t147 ^ 2;
t190 = t137 * t145;
t188 = t139 * t146;
t183 = qJD(1) * t146;
t181 = qJD(3) * t139;
t108 = t110 * t190 + 0.1e1;
t178 = 0.2e1 / t108 ^ 2 * (-t190 * t205 + (-t141 * t145 * t181 - t161) * t110);
t177 = 0.2e1 * t205;
t176 = 0.2e1 * t204;
t175 = -0.2e1 * t198;
t174 = t141 * t199;
t173 = t141 * t198;
t171 = t133 * t189;
t166 = t141 * t178;
t165 = 0.2e1 * t122 * t201;
t162 = t129 * t171;
t160 = t167 * t147;
t159 = t164 * t146;
t157 = t120 * t194 + t196;
t155 = t157 * t147;
t125 = -t138 * t188 + t186;
t124 = t138 * t147 + t139 * t185;
t114 = 0.1e1 / t117;
t106 = 0.1e1 / t108;
t102 = (t206 * t141 * t126 + t127 * t162) * t147;
t101 = t126 * t188 + t193 + (-t126 * t139 - t127 * t184) * t116;
t99 = t175 * t209 + (qJD(1) * t160 + 0.2e1 * t146 * t156) * t129;
t1 = [0.2e1 * t133 * t147 * t173 + (qJD(3) * t160 + t183 * t192) * t129, 0, t99, 0, 0, 0; (t109 * t166 + (t109 * t181 + (qJD(1) * t102 + t98) * t200) * t106) * t146 + (t110 * t166 * t102 + (-((-t103 * t162 - t206 * t181 - 0.2e1 * t173) * t126 + (t171 * t175 - t202 + (t202 + (-0.2e1 * t141 - t210) * t180) * t129) * t127) * t174 + (t110 * t181 + t141 * t177) * t102 + (-t109 + ((t144 - t145) * t137 * t133 * t129 * t127 + t206 * t172) * t110) * t141 * qJD(1)) * t106) * t147, 0 (t101 * t200 + t109 * t139) * t147 * t178 + ((t109 * t183 + (qJD(3) * t101 + t98) * t199) * t139 + (-t109 * t179 - (-t127 * t146 * t99 + t126 * t180 + t197 * t203 - t203 + (-qJD(3) * t126 - t127 * t182) * t116) * t174 + (t110 * t183 + t147 * t177) * t101 - ((-t99 + t182) * t126 + ((-0.1e1 + t197) * qJD(3) + (-t116 + t146) * t103) * t127) * t110 * t187) * t141) * t106, 0, 0, 0; (-t119 * t124 + t125 * t195) * t176 + (t125 * t165 - t119 * t138 * t159 + t207 * t196 + (t122 * t140 * t159 - t125 * t104 - t124 * t105 + t194 * t207) * t120) * t114, 0, t141 * t155 * t176 + (t155 * t181 + (t157 * t183 + ((qJD(6) * t119 + t165) * t138 + (-t104 * t138 + (-qJD(6) * t122 + t105) * t140) * t120) * t147) * t141) * t114, 0, 0, -0.2e1 * t204 + 0.2e1 * (t104 * t114 * t120 + (-t114 * t201 - t120 * t204) * t122) * t122;];
JaD_rot  = t1;
