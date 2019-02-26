% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:12
% EndTime: 2019-02-26 19:49:13
% DurationCPUTime: 0.82s
% Computational Cost: add. (2688->100), mult. (8196->220), div. (535->12), fcn. (10572->13), ass. (0->100)
t187 = sin(pkin(10));
t190 = cos(pkin(10));
t195 = cos(qJ(2));
t191 = cos(pkin(6));
t193 = sin(qJ(2));
t218 = t191 * t193;
t180 = t187 * t195 + t190 * t218;
t174 = t180 * qJD(2);
t188 = sin(pkin(6));
t194 = cos(qJ(4));
t192 = sin(qJ(4));
t217 = t191 * t195;
t205 = -t187 * t193 + t190 * t217;
t202 = t205 * t192;
t146 = qJD(4) * t202 + (qJD(4) * t188 * t190 + t174) * t194;
t222 = t188 * t192;
t165 = t190 * t222 - t194 * t205;
t162 = t165 ^ 2;
t219 = t188 * t195;
t183 = t191 * t192 + t194 * t219;
t178 = 0.1e1 / t183 ^ 2;
t157 = t162 * t178 + 0.1e1;
t155 = 0.1e1 / t157;
t184 = t191 * t194 - t192 * t219;
t221 = t188 * t193;
t208 = qJD(2) * t221;
t168 = qJD(4) * t184 - t194 * t208;
t177 = 0.1e1 / t183;
t227 = t165 * t178;
t127 = (t146 * t177 - t168 * t227) * t155;
t158 = atan2(t165, t183);
t153 = sin(t158);
t154 = cos(t158);
t207 = -t153 * t183 + t154 * t165;
t123 = t127 * t207 + t146 * t153 + t154 * t168;
t137 = t153 * t165 + t154 * t183;
t134 = 0.1e1 / t137;
t135 = 0.1e1 / t137 ^ 2;
t239 = t123 * t134 * t135;
t181 = t187 * t217 + t190 * t193;
t163 = -t181 * t194 + t187 * t222;
t238 = 0.2e1 * t163 * t239;
t225 = t168 * t177 * t178;
t237 = (t146 * t227 - t162 * t225) / t157 ^ 2;
t210 = t165 * t221;
t203 = t177 * t180 + t178 * t210;
t236 = t194 * t203;
t220 = t188 * t194;
t164 = t181 * t192 + t187 * t220;
t209 = t187 * t218;
t182 = t190 * t195 - t209;
t186 = sin(pkin(11));
t189 = cos(pkin(11));
t152 = t164 * t189 + t182 * t186;
t148 = 0.1e1 / t152;
t149 = 0.1e1 / t152 ^ 2;
t235 = t135 * t163;
t216 = qJD(2) * t195;
t176 = -qJD(2) * t209 + t190 * t216;
t143 = -qJD(4) * t163 + t176 * t192;
t175 = t181 * qJD(2);
t142 = t143 * t189 - t175 * t186;
t234 = t142 * t148 * t149;
t144 = qJD(4) * t164 - t176 * t194;
t233 = t144 * t135;
t232 = t148 * t186;
t151 = t164 * t186 - t182 * t189;
t231 = t149 * t151;
t230 = t151 * t189;
t229 = t153 * t163;
t228 = t154 * t163;
t226 = t165 * t184;
t224 = t182 * t192;
t223 = t182 * t194;
t215 = qJD(4) * t192;
t161 = t163 ^ 2;
t133 = t135 * t161 + 0.1e1;
t214 = 0.2e1 * (-t161 * t239 + t163 * t233) / t133 ^ 2;
t147 = t151 ^ 2;
t140 = t147 * t149 + 0.1e1;
t141 = t143 * t186 + t175 * t189;
t213 = 0.2e1 * (t141 * t231 - t147 * t234) / t140 ^ 2;
t211 = t151 * t234;
t166 = t190 * t220 + t202;
t206 = -t166 * t177 + t178 * t226;
t204 = qJD(4) * t223 - t175 * t192;
t173 = t205 * qJD(2);
t167 = -qJD(4) * t183 + t192 * t208;
t160 = -t181 * t186 + t189 * t224;
t159 = t181 * t189 + t186 * t224;
t145 = qJD(4) * t165 + t174 * t192;
t138 = 0.1e1 / t140;
t130 = 0.1e1 / t133;
t129 = t155 * t236;
t128 = t206 * t155;
t125 = (t153 * t180 - t154 * t221) * t194 + t207 * t129;
t124 = -t128 * t207 + t153 * t166 + t154 * t184;
t122 = 0.2e1 * t206 * t237 + (0.2e1 * t225 * t226 - t145 * t177 + (-t146 * t184 - t165 * t167 - t166 * t168) * t178) * t155;
t120 = -0.2e1 * t236 * t237 + (-t203 * t215 + (-0.2e1 * t210 * t225 + t173 * t177 + (-t168 * t180 + (t146 * t193 + t165 * t216) * t188) * t178) * t194) * t155;
t1 = [0, t120, 0, t122, 0, 0; 0 (t125 * t235 + t134 * t223) * t214 + ((t175 * t194 + t182 * t215) * t134 + (-t233 + t238) * t125 + (t223 * t123 - (t120 * t165 + t129 * t146 + (t193 * t215 - t194 * t216) * t188 + (-t129 * t183 + t180 * t194) * t127) * t228 - (-t180 * t215 - t120 * t183 - t129 * t168 + t173 * t194 + (-t129 * t165 + t193 * t220) * t127) * t229) * t135) * t130, 0 (t124 * t235 - t134 * t164) * t214 + (t124 * t238 + t143 * t134 + (-t164 * t123 - t124 * t144 - (t122 * t165 - t128 * t146 + t167 + (t128 * t183 + t166) * t127) * t228 - (-t122 * t183 + t128 * t168 - t145 + (t128 * t165 - t184) * t127) * t229) * t135) * t130, 0, 0; 0 (-t148 * t159 + t160 * t231) * t213 + ((t176 * t189 + t186 * t204) * t148 + 0.2e1 * t160 * t211 + (-t159 * t142 - (-t176 * t186 + t189 * t204) * t151 - t160 * t141) * t149) * t138, 0 (-t149 * t230 + t232) * t163 * t213 + (-0.2e1 * t163 * t189 * t211 - t144 * t232 + (t144 * t230 + (t141 * t189 + t142 * t186) * t163) * t149) * t138, 0, 0;];
JaD_rot  = t1;
