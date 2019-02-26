% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:02
% EndTime: 2019-02-26 21:02:03
% DurationCPUTime: 0.98s
% Computational Cost: add. (4034->125), mult. (6168->273), div. (1114->15), fcn. (7752->9), ass. (0->111)
t157 = cos(qJ(4));
t155 = sin(qJ(4));
t206 = qJ(1) + pkin(10);
t185 = sin(t206);
t180 = t185 * t155;
t147 = cos(t206);
t158 = cos(qJ(3));
t216 = t147 * t158;
t135 = t157 * t216 + t180;
t129 = 0.1e1 / t135 ^ 2;
t146 = t147 ^ 2;
t156 = sin(qJ(3));
t151 = t156 ^ 2;
t218 = t146 * t151;
t197 = t129 * t218;
t119 = 0.1e1 + t197;
t177 = qJD(1) * t185;
t173 = t158 * t177;
t176 = t185 * qJD(4);
t210 = qJD(3) * t156;
t114 = (-t173 + t176) * t157 + (-t157 * t210 + (-qJD(4) * t158 + qJD(1)) * t155) * t147;
t128 = 0.1e1 / t135;
t225 = t114 * t128 * t129;
t182 = t218 * t225;
t209 = qJD(3) * t158;
t231 = (-t182 + (t146 * t156 * t209 - t147 * t151 * t177) * t129) / t119 ^ 2;
t217 = t147 * t156;
t131 = t147 * t157 + t158 * t180;
t172 = t155 * t176;
t207 = qJD(4) * t157;
t189 = t147 * t207;
t193 = t147 * t210;
t113 = t131 * qJD(1) + t155 * t193 - t158 * t189 - t172;
t179 = t185 * t157;
t134 = t155 * t216 - t179;
t148 = 0.1e1 / t155;
t149 = 0.1e1 / t155 ^ 2;
t152 = 0.1e1 / t156;
t153 = 0.1e1 / t156 ^ 2;
t192 = t153 * t209;
t215 = t148 * t152;
t230 = (t149 * t152 * t207 + t148 * t192) * t134 + t113 * t215;
t212 = t156 * t155;
t124 = atan2(-t131, t212);
t121 = cos(t124);
t120 = sin(t124);
t224 = t120 * t131;
t112 = t121 * t212 - t224;
t109 = 0.1e1 / t112;
t110 = 0.1e1 / t112 ^ 2;
t229 = 0.2e1 * t134;
t126 = t131 ^ 2;
t214 = t149 * t153;
t125 = t126 * t214 + 0.1e1;
t122 = 0.1e1 / t125;
t170 = t155 * t209 + t156 * t207;
t195 = t131 * t214;
t181 = t156 * t185;
t174 = qJD(3) * t181;
t175 = t157 * t177;
t208 = qJD(4) * t155;
t211 = qJD(1) * t147;
t115 = -t155 * t174 - t147 * t208 - t175 + (t155 * t211 + t157 * t176) * t158;
t198 = t115 * t215;
t101 = (t170 * t195 - t198) * t122;
t166 = -t101 * t131 + t170;
t97 = (-t101 * t212 - t115) * t120 + t166 * t121;
t228 = t109 * t110 * t97;
t150 = t148 * t149;
t154 = t152 / t151;
t190 = t153 * t207;
t227 = (t115 * t195 + (-t149 * t154 * t209 - t150 * t190) * t126) / t125 ^ 2;
t226 = t110 * t134;
t223 = t120 * t134;
t222 = t120 * t156;
t221 = t121 * t131;
t220 = t121 * t134;
t219 = t121 * t158;
t213 = t149 * t157;
t127 = t134 ^ 2;
t107 = t110 * t127 + 0.1e1;
t205 = 0.2e1 / t107 ^ 2 * (-t113 * t226 - t127 * t228);
t204 = 0.2e1 * t228;
t203 = 0.2e1 * t231;
t202 = -0.2e1 * t227;
t201 = t152 * t227;
t200 = t110 * t223;
t196 = t131 * t215;
t194 = t148 * t153 * t158;
t191 = t157 * t209;
t188 = t109 * t205;
t187 = t110 * t205;
t186 = t134 * t204;
t184 = t217 * t229;
t183 = t148 * t201;
t169 = t131 * t194 + t185;
t108 = t169 * t122;
t178 = t185 - t108;
t133 = -t147 * t155 + t158 * t179;
t171 = t131 * t213 - t133 * t148;
t168 = t129 * t133 * t147 - t185 * t128;
t117 = 0.1e1 / t119;
t116 = t135 * qJD(1) - t157 * t174 - t158 * t172 - t189;
t105 = 0.1e1 / t107;
t104 = t171 * t152 * t122;
t100 = (-t120 + (t121 * t196 + t120) * t122) * t134;
t99 = -t108 * t221 + (t178 * t222 + t219) * t155;
t98 = t121 * t156 * t157 - t120 * t133 + (-t120 * t212 - t221) * t104;
t96 = t169 * t202 + (t115 * t194 + t211 + (-t149 * t158 * t190 + (-0.2e1 * t154 * t158 ^ 2 - t152) * t148 * qJD(3)) * t131) * t122;
t94 = -0.2e1 * t171 * t201 + (-t171 * t192 + (t115 * t213 - t116 * t148 + (t133 * t213 + (-0.2e1 * t150 * t157 ^ 2 - t148) * t131) * qJD(4)) * t152) * t122;
t1 = [t230 * t122 + t183 * t229, 0, t96, t94, 0, 0; t131 * t188 + (-t115 * t109 + (t100 * t113 + t131 * t97) * t110) * t105 + (t100 * t187 + (t100 * t204 + (t113 * t122 - t113 - (-t101 * t122 * t196 + t202) * t134) * t110 * t120 + (-(-0.2e1 * t131 * t183 - t101) * t226 + (-(t101 + t198) * t134 + t230 * t131) * t110 * t122) * t121) * t105) * t134, 0, t99 * t134 * t187 + (t99 * t186 + (-(-t96 * t221 + (t101 * t224 - t115 * t121) * t108) * t134 + t99 * t113) * t110 + (-t109 * t217 - (-t108 * t222 + t120 * t181 + t219) * t226) * t207) * t105 + (t188 * t217 + ((-t147 * qJD(3) * t109 - (t178 * qJD(3) - t101) * t200) * t158 + (t109 * t177 + (t147 * t97 - (-t96 + t211) * t223 - (t178 * t101 - qJD(3)) * t220) * t110) * t156) * t105) * t155 (-t109 * t135 + t98 * t226) * t205 + (t98 * t186 + t114 * t109 - (-t116 + (-t101 * t157 - t155 * t94) * t156 - t166 * t104) * t200 + (t98 * t113 - t135 * t97 - (t191 - t156 * t208 - t104 * t115 - t131 * t94 + (-t104 * t212 - t133) * t101) * t220) * t110) * t105, 0, 0; t168 * t156 * t203 + (-t168 * t209 + ((qJD(1) * t128 + 0.2e1 * t133 * t225) * t147 + (-t185 * t114 - t116 * t147 + t133 * t177) * t129) * t156) * t117, 0 (t128 * t216 + t157 * t197) * t203 + (0.2e1 * t157 * t182 + (t173 + t193) * t128 + ((t114 * t158 + 0.2e1 * t151 * t175) * t147 + (t151 * t208 - 0.2e1 * t156 * t191) * t146) * t129) * t117, t129 * t184 * t231 + (t184 * t225 + (t113 * t217 + (-t147 * t209 + t156 * t177) * t134) * t129) * t117, 0, 0;];
JaD_rot  = t1;
