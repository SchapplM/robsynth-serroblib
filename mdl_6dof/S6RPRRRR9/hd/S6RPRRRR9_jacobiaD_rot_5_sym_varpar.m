% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR9_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:14
% EndTime: 2019-02-26 21:19:15
% DurationCPUTime: 0.74s
% Computational Cost: add. (1381->94), mult. (2734->205), div. (498->12), fcn. (3199->9), ass. (0->97)
t150 = sin(qJ(3));
t143 = 0.1e1 / t150 ^ 2;
t152 = cos(qJ(3));
t147 = t152 ^ 2;
t197 = t143 * t147;
t153 = cos(qJ(1));
t215 = 0.2e1 * t153;
t214 = t152 * t197;
t190 = t153 * t152;
t134 = atan2(-t190, t150);
t132 = sin(t134);
t133 = cos(t134);
t121 = -t132 * t190 + t133 * t150;
t118 = 0.1e1 / t121;
t149 = qJ(4) + qJ(5);
t139 = sin(t149);
t140 = cos(t149);
t151 = sin(qJ(1));
t193 = t151 * t140;
t129 = t139 * t153 + t150 * t193;
t125 = 0.1e1 / t129;
t142 = 0.1e1 / t150;
t119 = 0.1e1 / t121 ^ 2;
t126 = 0.1e1 / t129 ^ 2;
t148 = t153 ^ 2;
t137 = t148 * t197 + 0.1e1;
t135 = 0.1e1 / t137;
t213 = t135 - 0.1e1;
t145 = t151 ^ 2;
t196 = t145 * t147;
t114 = t119 * t196 + 0.1e1;
t187 = qJD(1) * t153;
t168 = t147 * t151 * t187;
t185 = qJD(3) * t152;
t188 = qJD(1) * t152;
t177 = t151 * t188;
t184 = qJD(3) * t153;
t111 = ((t150 * t184 + t177) * t142 + t184 * t197) * t135;
t199 = t133 * t152;
t105 = (-t111 * t153 + qJD(3)) * t199 + (t177 + (-t111 + t184) * t150) * t132;
t210 = t105 * t118 * t119;
t212 = (-t196 * t210 + (-t145 * t150 * t185 + t168) * t119) / t114 ^ 2;
t141 = qJD(4) + qJD(5);
t171 = qJD(1) * t150 + t141;
t161 = t151 * t185 + t171 * t153;
t172 = t141 * t150 + qJD(1);
t165 = t140 * t172;
t109 = t161 * t139 + t151 * t165;
t191 = t153 * t140;
t194 = t151 * t139;
t128 = t150 * t194 - t191;
t124 = t128 ^ 2;
t117 = t124 * t126 + 0.1e1;
t201 = t126 * t128;
t166 = t139 * t172;
t110 = t161 * t140 - t151 * t166;
t209 = t110 * t125 * t126;
t211 = (t109 * t201 - t124 * t209) / t117 ^ 2;
t208 = t111 * t132;
t207 = t111 * t152;
t206 = t119 * t151;
t205 = t119 * t152;
t163 = qJD(3) * (-t152 - t214) * t142;
t204 = (-t143 * t168 + t148 * t163) / t137 ^ 2;
t175 = 0.1e1 + t197;
t123 = t175 * t153 * t135;
t203 = t123 * t153;
t202 = t125 * t139;
t200 = t128 * t140;
t198 = t142 * t147;
t195 = t150 * t153;
t192 = t151 * t152;
t189 = qJD(1) * t151;
t186 = qJD(3) * t150;
t183 = 0.2e1 * t211;
t182 = -0.2e1 * t210;
t181 = t152 * t212;
t180 = t119 * t192;
t179 = t152 * t204;
t178 = t135 * t198;
t176 = t152 * t187;
t174 = 0.2e1 * t128 * t209;
t173 = t204 * t215;
t170 = t153 * t178;
t169 = t213 * t152 * t132;
t167 = t175 * t151;
t164 = t126 * t200 - t202;
t162 = -t171 * t151 + t152 * t184;
t131 = t150 * t191 - t194;
t130 = t139 * t195 + t193;
t115 = 0.1e1 / t117;
t112 = 0.1e1 / t114;
t108 = (-t133 * t170 - t169) * t151;
t107 = t132 * t195 + t199 + (-t132 * t150 - t133 * t190) * t123;
t106 = -t175 * t173 + (-qJD(1) * t167 + t163 * t215) * t135;
t102 = -0.2e1 * t211 + 0.2e1 * (t109 * t115 * t126 + (-t115 * t209 - t126 * t211) * t128) * t128;
t1 = [-0.2e1 * t151 * t142 * t179 + (-qJD(3) * t167 + t142 * t176) * t135, 0, t106, 0, 0, 0; (0.2e1 * t118 * t181 + (t118 * t186 + (qJD(1) * t108 + t105) * t205) * t112) * t153 + (-0.2e1 * t119 * t181 * t108 + (((t111 * t170 + t213 * t186 + 0.2e1 * t179) * t132 + (t173 * t198 + t207 + (-t207 + (0.2e1 * t152 + t214) * t184) * t135) * t133) * t180 + (-t119 * t186 + t152 * t182) * t108 + (t118 + ((t145 - t148) * t133 * t178 - t153 * t169) * t119) * t188) * t112) * t151, 0, 0.2e1 * (-t107 * t205 - t118 * t150) * t151 * t212 + ((t118 * t187 + (-qJD(3) * t107 - t105) * t206) * t150 + (t151 * qJD(3) * t118 + (-t106 * t133 * t153 + t132 * t184 + t203 * t208 - t208 + (-qJD(3) * t132 + t133 * t189) * t123) * t180 + (t119 * t187 + t151 * t182) * t107 + ((-t106 - t189) * t132 + ((-0.1e1 + t203) * qJD(3) + (-t123 + t153) * t111) * t133) * t150 * t206) * t152) * t112, 0, 0, 0; (-t125 * t130 + t131 * t201) * t183 + (t131 * t174 + t153 * t125 * t165 + t162 * t202 + (t153 * t128 * t166 - t131 * t109 - t130 * t110 - t162 * t200) * t126) * t115, 0, t164 * t183 * t192 + (-t164 * t176 + (t164 * t186 + ((t125 * t141 + t174) * t140 + (-t109 * t140 + (t128 * t141 - t110) * t139) * t126) * t152) * t151) * t115, t102, t102, 0;];
JaD_rot  = t1;
