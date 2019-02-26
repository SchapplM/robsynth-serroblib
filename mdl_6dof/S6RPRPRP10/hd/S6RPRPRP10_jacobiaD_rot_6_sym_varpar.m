% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_jacobiaD_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:48:37
% EndTime: 2019-02-26 20:48:38
% DurationCPUTime: 1.14s
% Computational Cost: add. (2070->121), mult. (6168->269), div. (1114->15), fcn. (7752->9), ass. (0->110)
t150 = sin(qJ(5));
t152 = sin(qJ(1));
t153 = cos(qJ(5));
t154 = cos(qJ(3));
t155 = cos(qJ(1));
t205 = t154 * t155;
t131 = t150 * t152 - t153 * t205;
t151 = sin(qJ(3));
t208 = t151 * t153;
t119 = atan2(-t131, -t208);
t117 = sin(t119);
t118 = cos(t119);
t126 = t131 ^ 2;
t144 = 0.1e1 / t151 ^ 2;
t148 = 0.1e1 / t153 ^ 2;
t212 = t144 * t148;
t125 = t126 * t212 + 0.1e1;
t122 = 0.1e1 / t125;
t147 = 0.1e1 / t153;
t211 = t144 * t154;
t186 = t147 * t211;
t168 = t131 * t186 - t155;
t111 = t168 * t122;
t227 = t111 + t155;
t229 = t227 * t117 * t151 + t118 * t154;
t206 = t152 * t154;
t185 = t150 * t206;
t135 = t153 * t155 - t185;
t129 = 0.1e1 / t135 ^ 2;
t142 = t151 ^ 2;
t146 = t152 ^ 2;
t214 = t142 * t146;
t189 = t129 * t214;
t124 = 0.1e1 + t189;
t172 = qJD(1) * t154 + qJD(5);
t173 = qJD(5) * t154 + qJD(1);
t200 = qJD(3) * t152;
t180 = t151 * t200;
t207 = t152 * t153;
t114 = -t173 * t207 + (-t172 * t155 + t180) * t150;
t128 = 0.1e1 / t135;
t221 = t114 * t128 * t129;
t170 = t214 * t221;
t201 = qJD(3) * t151;
t181 = t146 * t201;
t202 = qJD(1) * t155;
t184 = t152 * t202;
t228 = (-t170 + (t142 * t184 + t154 * t181) * t129) / t124 ^ 2;
t209 = t151 * t152;
t183 = t154 * t202;
t197 = qJD(5) * t153;
t113 = -t153 * t183 - t155 * t197 + (t173 * t150 + t153 * t201) * t152;
t134 = t150 * t155 + t153 * t206;
t143 = 0.1e1 / t151;
t198 = qJD(5) * t150;
t178 = t148 * t198;
t199 = qJD(3) * t154;
t182 = t144 * t199;
t213 = t143 * t147;
t226 = -(-t143 * t178 + t147 * t182) * t134 - t113 * t213;
t220 = t117 * t131;
t112 = -t118 * t208 - t220;
t108 = 0.1e1 / t112;
t109 = 0.1e1 / t112 ^ 2;
t225 = 0.2e1 * t150;
t165 = t151 * t198 - t153 * t199;
t187 = t131 * t212;
t169 = t173 * t155;
t179 = t155 * t201;
t115 = t150 * t169 + (t172 * t152 + t179) * t153;
t190 = t115 * t213;
t101 = (t165 * t187 + t190) * t122;
t163 = t101 * t131 - t165;
t97 = (t101 * t208 - t115) * t117 - t163 * t118;
t224 = t108 * t109 * t97;
t145 = t143 / t142;
t149 = t147 * t148;
t223 = (t115 * t187 + (t144 * t149 * t198 - t145 * t148 * t199) * t126) / t125 ^ 2;
t222 = t109 * t134;
t219 = t117 * t134;
t217 = t118 * t131;
t216 = t118 * t134;
t210 = t148 * t150;
t203 = qJD(1) * t152;
t127 = t134 ^ 2;
t107 = t109 * t127 + 0.1e1;
t196 = 0.2e1 / t107 ^ 2 * (-t113 * t222 - t127 * t224);
t195 = 0.2e1 * t224;
t194 = 0.2e1 * t228;
t193 = t143 * t223;
t192 = t109 * t219;
t188 = t131 * t213;
t177 = t108 * t196;
t176 = t109 * t196;
t175 = t134 * t195;
t174 = 0.2e1 * t134 * t209;
t171 = t147 * t193;
t133 = t150 * t205 + t207;
t167 = t129 * t133 * t152 + t128 * t155;
t166 = t131 * t210 + t133 * t147;
t120 = 0.1e1 / t124;
t116 = -qJD(1) * t185 - t150 * t179 - t152 * t198 + t153 * t169;
t105 = 0.1e1 / t107;
t104 = t166 * t143 * t122;
t100 = (-t117 + (-t118 * t188 + t117) * t122) * t134;
t99 = t111 * t217 - t229 * t153;
t98 = t118 * t151 * t150 - t117 * t133 + (t117 * t208 - t217) * t104;
t96 = 0.2e1 * t168 * t223 + (-t115 * t186 - t203 + (-t178 * t211 + (0.2e1 * t145 * t154 ^ 2 + t143) * t147 * qJD(3)) * t131) * t122;
t94 = -0.2e1 * t166 * t193 + (-t166 * t182 + (t115 * t210 + t116 * t147 + (t133 * t210 + (0.2e1 * t149 * t150 ^ 2 + t147) * t131) * qJD(5)) * t143) * t122;
t1 = [t226 * t122 - 0.2e1 * t134 * t171, 0, t96, 0, t94, 0; t131 * t177 + (-t115 * t108 + (t100 * t113 + t131 * t97) * t109) * t105 + (t100 * t176 + (t100 * t195 + (t113 * t122 - t113 - (t101 * t122 * t188 - 0.2e1 * t223) * t134) * t109 * t117 + (-(0.2e1 * t131 * t171 - t101) * t222 + (-(t101 - t190) * t134 + t226 * t131) * t109 * t122) * t118) * t105) * t134, 0, t99 * t134 * t176 + (t99 * t175 + (-(-t96 * t217 - (t101 * t220 - t115 * t118) * t111) * t134 + t99 * t113) * t109 + (t108 * t209 - t229 * t222) * t198) * t105 + (t177 * t209 + ((-t108 * t200 - (-qJD(3) * t227 + t101) * t192) * t154 + (-t108 * t202 + (t152 * t97 - (t96 + t203) * t219 - (-t101 * t227 + qJD(3)) * t216) * t109) * t151) * t105) * t153, 0 (-t108 * t135 + t98 * t222) * t196 + (t98 * t175 + t114 * t108 - (-t116 + (-t101 * t150 + t153 * t94) * t151 + t163 * t104) * t192 + (t98 * t113 - t135 * t97 - (t150 * t199 + t151 * t197 - t104 * t115 - t131 * t94 + (t104 * t208 - t133) * t101) * t216) * t109) * t105, 0; t167 * t151 * t194 + (-t167 * t199 + ((qJD(1) * t128 + 0.2e1 * t133 * t221) * t152 + (-t116 * t152 + (-qJD(1) * t133 + t114) * t155) * t129) * t151) * t120, 0 (t128 * t206 - t150 * t189) * t194 + (-0.2e1 * t150 * t170 + (t180 - t183) * t128 + ((t114 * t152 + t181 * t225) * t154 + (t146 * t197 + t184 * t225) * t142) * t129) * t120, 0, t129 * t174 * t228 + (t174 * t221 + (t113 * t209 + (-t151 * t202 - t152 * t199) * t134) * t129) * t120, 0;];
JaD_rot  = t1;
