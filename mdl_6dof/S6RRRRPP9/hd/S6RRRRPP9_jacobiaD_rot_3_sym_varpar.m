% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP9_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobiaD_rot_3_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:06
% EndTime: 2019-02-26 22:30:07
% DurationCPUTime: 0.75s
% Computational Cost: add. (1479->91), mult. (4303->201), div. (668->14), fcn. (5516->11), ass. (0->91)
t171 = sin(qJ(2));
t172 = sin(qJ(1));
t225 = cos(pkin(6));
t195 = t172 * t225;
t193 = t171 * t195;
t174 = cos(qJ(2));
t175 = cos(qJ(1));
t209 = t175 * t174;
t157 = -t193 + t209;
t170 = sin(qJ(3));
t173 = cos(qJ(3));
t169 = sin(pkin(6));
t213 = t169 * t172;
t185 = -t157 * t170 + t173 * t213;
t227 = t185 * qJD(3);
t194 = t175 * t225;
t192 = t174 * t194;
t210 = t172 * t171;
t153 = -t192 + t210;
t212 = t169 * t174;
t147 = atan2(-t153, -t212);
t145 = sin(t147);
t146 = cos(t147);
t151 = t153 ^ 2;
t165 = 0.1e1 / t169 ^ 2;
t167 = 0.1e1 / t174 ^ 2;
t150 = t151 * t165 * t167 + 0.1e1;
t148 = 0.1e1 / t150;
t164 = 0.1e1 / t169;
t166 = 0.1e1 / t174;
t199 = t153 * t164 * t166;
t226 = (t146 * t199 - t145) * t148 + t145;
t129 = -t145 * t153 - t146 * t212;
t126 = 0.1e1 / t129;
t144 = t157 * t173 + t170 * t213;
t138 = 0.1e1 / t144;
t127 = 0.1e1 / t129 ^ 2;
t139 = 0.1e1 / t144 ^ 2;
t182 = -t171 * t194 - t172 * t174;
t183 = -t175 * t171 - t174 * t195;
t135 = -t183 * qJD(1) - t182 * qJD(2);
t207 = qJD(2) * t171;
t196 = t167 * t207;
t184 = t135 * t166 + t153 * t196;
t215 = t148 * t164;
t118 = t184 * t215;
t188 = t145 * t212 - t146 * t153;
t200 = t146 * t169 * t171;
t114 = qJD(2) * t200 + t188 * t118 - t145 * t135;
t224 = t114 * t126 * t127;
t214 = t167 * t171;
t187 = t153 * t214 - t166 * t182;
t119 = t187 * t215;
t115 = t188 * t119 + t145 * t182 + t200;
t223 = t115 * t183;
t134 = t182 * qJD(1) + t183 * qJD(2);
t208 = qJD(1) * t169;
t197 = t175 * t208;
t124 = t144 * qJD(3) + t134 * t170 - t173 * t197;
t137 = t185 ^ 2;
t132 = t137 * t139 + 0.1e1;
t218 = t139 * t185;
t125 = t134 * t173 + t170 * t197 + t227;
t220 = t125 * t138 * t139;
t222 = (-t124 * t218 - t137 * t220) / t132 ^ 2;
t168 = t166 * t167;
t221 = (t135 * t153 * t167 + t151 * t168 * t207) * t165 / t150 ^ 2;
t191 = qJD(2) * t225 + qJD(1);
t206 = qJD(2) * t174;
t133 = -qJD(1) * t192 - t175 * t206 + t191 * t210;
t219 = t133 * t127;
t217 = t145 * t183;
t216 = t146 * t183;
t211 = t169 * t175;
t152 = t183 ^ 2;
t122 = t152 * t127 + 0.1e1;
t205 = 0.2e1 * (-t152 * t224 + t183 * t219) / t122 ^ 2;
t204 = 0.2e1 * t224;
t203 = 0.2e1 * t222;
t202 = -0.2e1 * t221;
t201 = t185 * t220;
t198 = t172 * t208;
t189 = t170 * t138 + t173 * t218;
t186 = -t170 * t182 + t173 * t211;
t142 = t170 * t211 + t173 * t182;
t136 = -qJD(1) * t193 - t172 * t207 + t191 * t209;
t130 = 0.1e1 / t132;
t120 = 0.1e1 / t122;
t117 = t226 * t183;
t113 = (t187 * t202 + (t135 * t214 + t136 * t166 + (-t182 * t214 + (0.2e1 * t168 * t171 ^ 2 + t166) * t153) * qJD(2)) * t148) * t164;
t1 = [(-t183 * t166 * t202 + (-t133 * t166 - t183 * t196) * t148) * t164, t113, 0, 0, 0, 0; t153 * t126 * t205 + (-t135 * t126 + (t114 * t153 + t117 * t133) * t127) * t120 - ((t117 * t204 - t226 * t219) * t120 + (t117 * t205 + ((t118 * t148 * t199 + t202) * t217 + (0.2e1 * t199 * t221 - t118 + (-t184 * t164 + t118) * t148) * t216) * t120) * t127) * t183 (-t126 * t157 - t127 * t223) * t205 + (-t204 * t223 + t134 * t126 + (-t157 * t114 + t115 * t133 + (t169 * t206 - t113 * t153 - t119 * t135 + (t119 * t212 + t182) * t118) * t216 + (t118 * t119 * t153 - t136 + (t113 * t174 + (-qJD(2) * t119 - t118) * t171) * t169) * t217) * t127) * t120, 0, 0, 0, 0; (t138 * t186 - t142 * t218) * t203 + ((t142 * qJD(3) - t136 * t170 + t173 * t198) * t138 - 0.2e1 * t142 * t201 + (t186 * t125 + (t186 * qJD(3) - t136 * t173 - t170 * t198) * t185 - t142 * t124) * t139) * t130, -t189 * t183 * t203 + (t189 * t133 - ((-qJD(3) * t138 + 0.2e1 * t201) * t173 + (t124 * t173 + (t125 + t227) * t170) * t139) * t183) * t130, -0.2e1 * t222 - 0.2e1 * (t124 * t139 * t130 - (-t130 * t220 - t139 * t222) * t185) * t185, 0, 0, 0;];
JaD_rot  = t1;
