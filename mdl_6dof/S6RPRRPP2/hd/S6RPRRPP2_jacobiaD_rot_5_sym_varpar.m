% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:54
% EndTime: 2019-02-26 20:56:55
% DurationCPUTime: 1.03s
% Computational Cost: add. (4034->125), mult. (6168->274), div. (1114->15), fcn. (7752->9), ass. (0->110)
t155 = sin(qJ(4));
t157 = cos(qJ(4));
t205 = qJ(1) + pkin(9);
t184 = sin(t205);
t147 = cos(t205);
t158 = cos(qJ(3));
t215 = t147 * t158;
t135 = t184 * t155 + t157 * t215;
t129 = 0.1e1 / t135 ^ 2;
t146 = t147 ^ 2;
t156 = sin(qJ(3));
t151 = t156 ^ 2;
t217 = t146 * t151;
t196 = t129 * t217;
t119 = 0.1e1 + t196;
t177 = qJD(1) * t184;
t173 = t158 * t177;
t176 = t184 * qJD(4);
t209 = qJD(3) * t156;
t114 = (-t173 + t176) * t157 + (-t157 * t209 + (-qJD(4) * t158 + qJD(1)) * t155) * t147;
t128 = 0.1e1 / t135;
t224 = t114 * t128 * t129;
t181 = t217 * t224;
t208 = qJD(3) * t158;
t230 = (-t181 + (t146 * t156 * t208 - t147 * t151 * t177) * t129) / t119 ^ 2;
t216 = t147 * t156;
t179 = t184 * t158;
t131 = t147 * t157 + t155 * t179;
t172 = t155 * t176;
t206 = qJD(4) * t157;
t188 = t147 * t206;
t192 = t147 * t209;
t113 = t131 * qJD(1) + t155 * t192 - t158 * t188 - t172;
t134 = t155 * t215 - t184 * t157;
t148 = 0.1e1 / t155;
t149 = 0.1e1 / t155 ^ 2;
t152 = 0.1e1 / t156;
t153 = 0.1e1 / t156 ^ 2;
t191 = t153 * t208;
t214 = t148 * t152;
t229 = (t149 * t152 * t206 + t148 * t191) * t134 + t113 * t214;
t211 = t156 * t155;
t124 = atan2(-t131, t211);
t121 = cos(t124);
t120 = sin(t124);
t223 = t120 * t131;
t112 = t121 * t211 - t223;
t109 = 0.1e1 / t112;
t110 = 0.1e1 / t112 ^ 2;
t228 = 0.2e1 * t134;
t126 = t131 ^ 2;
t213 = t149 * t153;
t125 = t126 * t213 + 0.1e1;
t122 = 0.1e1 / t125;
t170 = t155 * t208 + t156 * t206;
t194 = t131 * t213;
t180 = t184 * t156;
t174 = qJD(3) * t180;
t175 = t157 * t177;
t207 = qJD(4) * t155;
t210 = qJD(1) * t147;
t115 = -t155 * t174 - t147 * t207 - t175 + (t155 * t210 + t157 * t176) * t158;
t197 = t115 * t214;
t101 = (t170 * t194 - t197) * t122;
t166 = -t101 * t131 + t170;
t97 = (-t101 * t211 - t115) * t120 + t166 * t121;
t227 = t109 * t110 * t97;
t150 = t148 * t149;
t154 = t152 / t151;
t189 = t153 * t206;
t226 = (t115 * t194 + (-t149 * t154 * t208 - t150 * t189) * t126) / t125 ^ 2;
t225 = t110 * t134;
t222 = t120 * t134;
t221 = t120 * t156;
t220 = t121 * t131;
t219 = t121 * t134;
t218 = t121 * t158;
t212 = t149 * t157;
t127 = t134 ^ 2;
t107 = t110 * t127 + 0.1e1;
t204 = 0.2e1 / t107 ^ 2 * (-t113 * t225 - t127 * t227);
t203 = 0.2e1 * t227;
t202 = 0.2e1 * t230;
t201 = -0.2e1 * t226;
t200 = t152 * t226;
t199 = t110 * t222;
t195 = t131 * t214;
t193 = t148 * t153 * t158;
t190 = t157 * t208;
t187 = t109 * t204;
t186 = t110 * t204;
t185 = t134 * t203;
t183 = t216 * t228;
t182 = t148 * t200;
t169 = t131 * t193 + t184;
t108 = t169 * t122;
t178 = t184 - t108;
t133 = -t147 * t155 + t157 * t179;
t171 = t131 * t212 - t133 * t148;
t168 = t129 * t133 * t147 - t184 * t128;
t117 = 0.1e1 / t119;
t116 = t135 * qJD(1) - t157 * t174 - t158 * t172 - t188;
t105 = 0.1e1 / t107;
t104 = t171 * t152 * t122;
t100 = (-t120 + (t121 * t195 + t120) * t122) * t134;
t99 = -t108 * t220 + (t178 * t221 + t218) * t155;
t98 = t121 * t156 * t157 - t120 * t133 + (-t120 * t211 - t220) * t104;
t96 = t169 * t201 + (t115 * t193 + t210 + (-t149 * t158 * t189 + (-0.2e1 * t154 * t158 ^ 2 - t152) * t148 * qJD(3)) * t131) * t122;
t94 = -0.2e1 * t171 * t200 + (-t171 * t191 + (t115 * t212 - t116 * t148 + (t133 * t212 + (-0.2e1 * t150 * t157 ^ 2 - t148) * t131) * qJD(4)) * t152) * t122;
t1 = [t229 * t122 + t182 * t228, 0, t96, t94, 0, 0; t131 * t187 + (-t115 * t109 + (t100 * t113 + t131 * t97) * t110) * t105 + (t100 * t186 + (t100 * t203 + (t113 * t122 - t113 - (-t101 * t122 * t195 + t201) * t134) * t110 * t120 + (-(-0.2e1 * t131 * t182 - t101) * t225 + (-(t101 + t197) * t134 + t229 * t131) * t110 * t122) * t121) * t105) * t134, 0, t99 * t134 * t186 + (t99 * t185 + (-(-t96 * t220 + (t101 * t223 - t115 * t121) * t108) * t134 + t99 * t113) * t110 + (-t109 * t216 - (-t108 * t221 + t120 * t180 + t218) * t225) * t206) * t105 + (t187 * t216 + ((-t147 * qJD(3) * t109 - (t178 * qJD(3) - t101) * t199) * t158 + (t109 * t177 + (t147 * t97 - (-t96 + t210) * t222 - (t178 * t101 - qJD(3)) * t219) * t110) * t156) * t105) * t155 (-t109 * t135 + t98 * t225) * t204 + (t98 * t185 + t114 * t109 - (-t116 + (-t101 * t157 - t155 * t94) * t156 - t166 * t104) * t199 + (t98 * t113 - t135 * t97 - (t190 - t156 * t207 - t104 * t115 - t131 * t94 + (-t104 * t211 - t133) * t101) * t219) * t110) * t105, 0, 0; t168 * t156 * t202 + (-t168 * t208 + ((qJD(1) * t128 + 0.2e1 * t133 * t224) * t147 + (-t184 * t114 - t116 * t147 + t133 * t177) * t129) * t156) * t117, 0 (t128 * t215 + t157 * t196) * t202 + (0.2e1 * t157 * t181 + (t173 + t192) * t128 + ((t114 * t158 + 0.2e1 * t151 * t175) * t147 + (t151 * t207 - 0.2e1 * t156 * t190) * t146) * t129) * t117, t129 * t183 * t230 + (t183 * t224 + (t113 * t216 + (-t147 * t208 + t156 * t177) * t134) * t129) * t117, 0, 0;];
JaD_rot  = t1;
