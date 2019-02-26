% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPP2_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_jacobiaD_rot_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:25:31
% EndTime: 2019-02-26 22:25:32
% DurationCPUTime: 0.79s
% Computational Cost: add. (3645->97), mult. (3810->208), div. (753->12), fcn. (4455->9), ass. (0->95)
t173 = sin(qJ(1));
t235 = 0.2e1 * t173;
t169 = t173 ^ 2;
t171 = qJ(2) + qJ(3);
t166 = sin(t171);
t162 = t166 ^ 2;
t167 = cos(t171);
t164 = 0.1e1 / t167 ^ 2;
t220 = t162 * t164;
t157 = t169 * t220 + 0.1e1;
t155 = 0.1e1 / t157;
t163 = 0.1e1 / t167;
t175 = cos(qJ(1));
t206 = qJD(1) * t175;
t196 = t166 * t206;
t168 = qJD(2) + qJD(3);
t214 = t168 * t173;
t199 = t164 * t214;
t129 = (-(-t167 * t214 - t196) * t163 + t162 * t199) * t155;
t234 = t129 - t214;
t174 = cos(qJ(4));
t208 = t174 * t175;
t172 = sin(qJ(4));
t210 = t173 * t172;
t151 = t167 * t208 + t210;
t211 = t173 * t166;
t154 = atan2(-t211, -t167);
t153 = cos(t154);
t152 = sin(t154);
t200 = t152 * t211;
t139 = -t153 * t167 - t200;
t136 = 0.1e1 / t139;
t145 = 0.1e1 / t151;
t137 = 0.1e1 / t139 ^ 2;
t146 = 0.1e1 / t151 ^ 2;
t233 = t155 - 0.1e1;
t222 = t153 * t166;
t124 = (-t129 * t173 + t168) * t222 + (t234 * t167 - t196) * t152;
t232 = t124 * t136 * t137;
t184 = t167 * t210 + t208;
t213 = t168 * t175;
t197 = t166 * t213;
t133 = t184 * qJD(1) - t151 * qJD(4) + t172 * t197;
t209 = t173 * t174;
t212 = t172 * t175;
t150 = t167 * t212 - t209;
t144 = t150 ^ 2;
t143 = t144 * t146 + 0.1e1;
t225 = t146 * t150;
t190 = -qJD(1) * t167 + qJD(4);
t191 = qJD(4) * t167 - qJD(1);
t134 = -t191 * t212 + (t173 * t190 - t197) * t174;
t230 = t134 * t145 * t146;
t231 = (-t133 * t225 - t144 * t230) / t143 ^ 2;
t161 = t166 * t162;
t217 = t163 * t166;
t183 = t168 * (t161 * t163 * t164 + t217);
t218 = t162 * t173;
t188 = t206 * t218;
t229 = (t164 * t188 + t169 * t183) / t157 ^ 2;
t228 = t137 * t166;
t227 = t137 * t175;
t226 = t145 * t172;
t224 = t150 * t174;
t223 = t152 * t173;
t221 = t162 * t163;
t170 = t175 ^ 2;
t219 = t162 * t170;
t216 = t166 * t175;
t215 = t167 * t168;
t207 = qJD(1) * t173;
t132 = t137 * t219 + 0.1e1;
t205 = 0.2e1 * (-t219 * t232 + (t166 * t170 * t215 - t188) * t137) / t132 ^ 2;
t204 = 0.2e1 * t232;
t203 = -0.2e1 * t231;
t202 = t137 * t216;
t201 = t150 * t230;
t195 = 0.1e1 + t220;
t194 = t166 * t205;
t193 = -0.2e1 * t166 * t229;
t192 = t229 * t235;
t189 = t153 * t155 * t221;
t187 = t195 * t175;
t186 = t190 * t175;
t185 = t146 * t224 - t226;
t149 = -t167 * t209 + t212;
t141 = 0.1e1 / t143;
t140 = t195 * t173 * t155;
t130 = 0.1e1 / t132;
t128 = (t152 * t166 * t233 - t173 * t189) * t175;
t126 = -t167 * t223 + t222 + (t152 * t167 - t153 * t211) * t140;
t125 = -t195 * t192 + (qJD(1) * t187 + t183 * t235) * t155;
t122 = t185 * t203 * t216 + (t185 * t167 * t213 + (-t185 * t207 + ((-qJD(4) * t145 - 0.2e1 * t201) * t174 + (-t133 * t174 + (-qJD(4) * t150 + t134) * t172) * t146) * t175) * t166) * t141;
t121 = (t126 * t228 - t136 * t167) * t175 * t205 + ((-t136 * t207 + (-t126 * t168 - t124) * t227) * t167 + (-t136 * t213 - (-t125 * t153 * t173 - t234 * t152 + (t129 * t223 - t152 * t168 - t153 * t206) * t140) * t202 + (t137 * t207 + t175 * t204) * t126 - ((t125 - t206) * t152 + ((-t140 * t173 + 0.1e1) * t168 + (t140 - t173) * t129) * t153) * t167 * t227) * t166) * t130;
t1 = [t163 * t175 * t193 + (t168 * t187 - t207 * t217) * t155, t125, t125, 0, 0, 0; (t136 * t194 + (-t136 * t215 + (qJD(1) * t128 + t124) * t228) * t130) * t173 + (t137 * t194 * t128 + (-((t193 - t215 + (t129 * t163 * t218 + t215) * t155) * t152 + (t192 * t221 - t129 * t166 + (-t161 * t199 + (t129 - 0.2e1 * t214) * t166) * t155) * t153) * t202 + (-t137 * t215 + t166 * t204) * t128 + (-t136 + ((-t169 + t170) * t189 + t233 * t200) * t137) * t166 * qJD(1)) * t130) * t175, t121, t121, 0, 0, 0; 0.2e1 * (t145 * t184 + t149 * t225) * t231 + (0.2e1 * t149 * t201 - t191 * t145 * t209 + (t168 * t211 + t186) * t226 + (t149 * t133 + t184 * t134 - t186 * t224 - (t166 * t168 * t174 + t172 * t191) * t150 * t173) * t146) * t141, t122, t122, t203 + 0.2e1 * (-t133 * t141 * t146 + (-t141 * t230 - t146 * t231) * t150) * t150, 0, 0;];
JaD_rot  = t1;
