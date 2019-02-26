% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:48
% EndTime: 2019-02-26 21:30:49
% DurationCPUTime: 0.59s
% Computational Cost: add. (1334->89), mult. (4303->202), div. (668->14), fcn. (5516->11), ass. (0->93)
t171 = sin(qJ(2));
t172 = sin(qJ(1));
t174 = cos(qJ(2));
t175 = cos(qJ(1));
t222 = cos(pkin(6));
t191 = t175 * t222;
t154 = t171 * t191 + t172 * t174;
t192 = t172 * t222;
t155 = t175 * t171 + t174 * t192;
t135 = t155 * qJD(1) + t154 * qJD(2);
t188 = t174 * t191;
t207 = t172 * t171;
t153 = -t188 + t207;
t151 = t153 ^ 2;
t169 = sin(pkin(6));
t165 = 0.1e1 / t169 ^ 2;
t167 = 0.1e1 / t174 ^ 2;
t150 = t151 * t165 * t167 + 0.1e1;
t166 = 0.1e1 / t174;
t168 = t166 * t167;
t204 = qJD(2) * t171;
t217 = (t135 * t153 * t167 + t151 * t168 * t204) * t165 / t150 ^ 2;
t225 = -0.2e1 * t217;
t164 = 0.1e1 / t169;
t224 = t153 * t164;
t209 = t169 * t174;
t149 = atan2(t153, t209);
t145 = sin(t149);
t146 = cos(t149);
t147 = 0.1e1 / t150;
t196 = t166 * t224;
t223 = (t146 * t196 - t145) * t147 + t145;
t129 = t145 * t153 + t146 * t209;
t126 = 0.1e1 / t129;
t189 = t171 * t192;
t206 = t175 * t174;
t157 = -t189 + t206;
t170 = sin(qJ(5));
t173 = cos(qJ(5));
t210 = t169 * t172;
t144 = t157 * t173 - t170 * t210;
t138 = 0.1e1 / t144;
t127 = 0.1e1 / t129 ^ 2;
t139 = 0.1e1 / t144 ^ 2;
t152 = t155 ^ 2;
t122 = t152 * t127 + 0.1e1;
t187 = qJD(2) * t222 + qJD(1);
t203 = qJD(2) * t174;
t133 = -qJD(1) * t188 - t175 * t203 + t187 * t207;
t215 = t133 * t127;
t193 = t167 * t204;
t182 = (t135 * t166 + t153 * t193) * t164;
t118 = t147 * t182;
t184 = -t145 * t209 + t146 * t153;
t197 = t146 * t169 * t171;
t114 = -qJD(2) * t197 + t184 * t118 + t145 * t135;
t220 = t114 * t126 * t127;
t221 = (-t152 * t220 - t155 * t215) / t122 ^ 2;
t211 = t167 * t171;
t183 = t153 * t211 + t154 * t166;
t119 = t183 * t164 * t147;
t115 = t184 * t119 + t145 * t154 - t197;
t219 = t115 * t155;
t134 = t154 * qJD(1) + t155 * qJD(2);
t205 = qJD(1) * t169;
t194 = t175 * t205;
t124 = t144 * qJD(5) - t134 * t170 + t173 * t194;
t143 = t157 * t170 + t173 * t210;
t137 = t143 ^ 2;
t132 = t137 * t139 + 0.1e1;
t214 = t139 * t143;
t202 = qJD(5) * t143;
t125 = -t134 * t173 - t170 * t194 - t202;
t216 = t125 * t138 * t139;
t218 = (t124 * t214 - t137 * t216) / t132 ^ 2;
t213 = t145 * t155;
t212 = t146 * t155;
t208 = t169 * t175;
t201 = -0.2e1 * t221;
t200 = -0.2e1 * t220;
t199 = 0.2e1 * t218;
t198 = t143 * t216;
t195 = t172 * t205;
t190 = t166 * t225;
t185 = t170 * t138 - t173 * t214;
t141 = -t154 * t170 + t173 * t208;
t142 = -t154 * t173 - t170 * t208;
t136 = -qJD(1) * t189 - t172 * t204 + t187 * t206;
t130 = 0.1e1 / t132;
t120 = 0.1e1 / t122;
t117 = t223 * t155;
t113 = (t183 * t225 + (t135 * t211 + t136 * t166 + (t154 * t211 + (0.2e1 * t168 * t171 ^ 2 + t166) * t153) * qJD(2)) * t147) * t164;
t1 = [(t155 * t190 + (-t133 * t166 + t155 * t193) * t147) * t164, t113, 0, 0, 0, 0; t153 * t126 * t201 + (t135 * t126 + (-t114 * t153 - t117 * t133) * t127) * t120 + ((t117 * t200 - t223 * t215) * t120 + (t117 * t201 + ((-t118 * t147 * t196 + 0.2e1 * t217) * t213 + (t190 * t224 + t118 + (-t118 + t182) * t147) * t212) * t120) * t127) * t155, 0.2e1 * (t126 * t157 - t127 * t219) * t221 + (t200 * t219 + t134 * t126 + (t157 * t114 - t115 * t133 + (-t169 * t203 + t113 * t153 + t119 * t135 + (-t119 * t209 + t154) * t118) * t212 + (-t118 * t119 * t153 + t136 + (-t113 * t174 + (qJD(2) * t119 + t118) * t171) * t169) * t213) * t127) * t120, 0, 0, 0, 0; (-t138 * t141 + t142 * t214) * t199 + ((t142 * qJD(5) - t136 * t170 - t173 * t195) * t138 + 0.2e1 * t142 * t198 + (-t141 * t125 - (-t141 * qJD(5) - t136 * t173 + t170 * t195) * t143 - t142 * t124) * t139) * t130, t185 * t155 * t199 + (t185 * t133 + ((-qJD(5) * t138 - 0.2e1 * t198) * t173 + (t124 * t173 + (t125 - t202) * t170) * t139) * t155) * t130, 0, 0, -0.2e1 * t218 + 0.2e1 * (t124 * t139 * t130 + (-t130 * t216 - t139 * t218) * t143) * t143, 0;];
JaD_rot  = t1;
