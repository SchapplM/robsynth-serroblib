% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:59
% EndTime: 2019-02-26 20:54:00
% DurationCPUTime: 0.72s
% Computational Cost: add. (1820->94), mult. (2734->205), div. (498->12), fcn. (3199->9), ass. (0->97)
t153 = sin(qJ(3));
t147 = 0.1e1 / t153 ^ 2;
t155 = cos(qJ(3));
t151 = t155 ^ 2;
t200 = t147 * t151;
t156 = cos(qJ(1));
t218 = 0.2e1 * t156;
t217 = t155 * t200;
t193 = t156 * t155;
t137 = atan2(-t193, t153);
t135 = sin(t137);
t136 = cos(t137);
t124 = -t135 * t193 + t136 * t153;
t121 = 0.1e1 / t124;
t144 = pkin(10) + qJ(5) + qJ(6);
t142 = sin(t144);
t143 = cos(t144);
t154 = sin(qJ(1));
t196 = t154 * t143;
t132 = t142 * t156 + t153 * t196;
t128 = 0.1e1 / t132;
t146 = 0.1e1 / t153;
t122 = 0.1e1 / t124 ^ 2;
t129 = 0.1e1 / t132 ^ 2;
t152 = t156 ^ 2;
t140 = t152 * t200 + 0.1e1;
t138 = 0.1e1 / t140;
t216 = t138 - 0.1e1;
t149 = t154 ^ 2;
t199 = t149 * t151;
t120 = t122 * t199 + 0.1e1;
t190 = qJD(1) * t156;
t171 = t151 * t154 * t190;
t188 = qJD(3) * t155;
t191 = qJD(1) * t155;
t180 = t154 * t191;
t187 = qJD(3) * t156;
t114 = ((t153 * t187 + t180) * t146 + t187 * t200) * t138;
t202 = t136 * t155;
t108 = (-t114 * t156 + qJD(3)) * t202 + (t180 + (-t114 + t187) * t153) * t135;
t213 = t108 * t121 * t122;
t215 = (-t199 * t213 + (-t149 * t153 * t188 + t171) * t122) / t120 ^ 2;
t145 = qJD(5) + qJD(6);
t174 = qJD(1) * t153 + t145;
t164 = t154 * t188 + t174 * t156;
t175 = t145 * t153 + qJD(1);
t168 = t143 * t175;
t112 = t164 * t142 + t154 * t168;
t194 = t156 * t143;
t197 = t154 * t142;
t131 = t153 * t197 - t194;
t127 = t131 ^ 2;
t117 = t127 * t129 + 0.1e1;
t204 = t129 * t131;
t169 = t142 * t175;
t113 = t164 * t143 - t154 * t169;
t212 = t113 * t128 * t129;
t214 = (t112 * t204 - t127 * t212) / t117 ^ 2;
t211 = t114 * t135;
t210 = t114 * t155;
t209 = t122 * t154;
t208 = t122 * t155;
t166 = qJD(3) * (-t155 - t217) * t146;
t207 = (-t147 * t171 + t152 * t166) / t140 ^ 2;
t178 = 0.1e1 + t200;
t126 = t178 * t156 * t138;
t206 = t126 * t156;
t205 = t128 * t142;
t203 = t131 * t143;
t201 = t146 * t151;
t198 = t153 * t156;
t195 = t154 * t155;
t192 = qJD(1) * t154;
t189 = qJD(3) * t153;
t186 = 0.2e1 * t214;
t185 = -0.2e1 * t213;
t184 = t155 * t215;
t183 = t122 * t195;
t182 = t155 * t207;
t181 = t138 * t201;
t179 = t155 * t190;
t177 = 0.2e1 * t131 * t212;
t176 = t207 * t218;
t173 = t156 * t181;
t172 = t216 * t155 * t135;
t170 = t178 * t154;
t167 = t129 * t203 - t205;
t165 = -t174 * t154 + t155 * t187;
t134 = t153 * t194 - t197;
t133 = t142 * t198 + t196;
t118 = 0.1e1 / t120;
t115 = 0.1e1 / t117;
t111 = (-t136 * t173 - t172) * t154;
t110 = t135 * t198 + t202 + (-t135 * t153 - t136 * t193) * t126;
t109 = -t178 * t176 + (-qJD(1) * t170 + t166 * t218) * t138;
t105 = -0.2e1 * t214 + 0.2e1 * (t112 * t115 * t129 + (-t115 * t212 - t129 * t214) * t131) * t131;
t1 = [-0.2e1 * t154 * t146 * t182 + (-qJD(3) * t170 + t146 * t179) * t138, 0, t109, 0, 0, 0; (0.2e1 * t121 * t184 + (t121 * t189 + (qJD(1) * t111 + t108) * t208) * t118) * t156 + (-0.2e1 * t122 * t184 * t111 + (((t114 * t173 + t216 * t189 + 0.2e1 * t182) * t135 + (t176 * t201 + t210 + (-t210 + (0.2e1 * t155 + t217) * t187) * t138) * t136) * t183 + (-t122 * t189 + t155 * t185) * t111 + (t121 + ((t149 - t152) * t136 * t181 - t156 * t172) * t122) * t191) * t118) * t154, 0, 0.2e1 * (-t110 * t208 - t121 * t153) * t154 * t215 + ((t121 * t190 + (-qJD(3) * t110 - t108) * t209) * t153 + (t154 * qJD(3) * t121 + (-t109 * t136 * t156 + t135 * t187 + t206 * t211 - t211 + (-qJD(3) * t135 + t136 * t192) * t126) * t183 + (t122 * t190 + t154 * t185) * t110 + ((-t109 - t192) * t135 + ((-0.1e1 + t206) * qJD(3) + (-t126 + t156) * t114) * t136) * t153 * t209) * t155) * t118, 0, 0, 0; (-t128 * t133 + t134 * t204) * t186 + (t134 * t177 + t156 * t128 * t168 + t165 * t205 + (t156 * t131 * t169 - t134 * t112 - t133 * t113 - t165 * t203) * t129) * t115, 0, t167 * t186 * t195 + (-t167 * t179 + (t167 * t189 + ((t128 * t145 + t177) * t143 + (-t112 * t143 + (t131 * t145 - t113) * t142) * t129) * t155) * t154) * t115, 0, t105, t105;];
JaD_rot  = t1;
