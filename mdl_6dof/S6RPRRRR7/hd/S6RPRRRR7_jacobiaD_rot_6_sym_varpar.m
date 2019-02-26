% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR7
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
% Datum: 2019-02-26 21:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:18:11
% EndTime: 2019-02-26 21:18:12
% DurationCPUTime: 0.76s
% Computational Cost: add. (7389->95), mult. (5101->207), div. (1026->12), fcn. (5942->9), ass. (0->97)
t174 = cos(qJ(1));
t234 = 0.2e1 * t174;
t170 = t174 ^ 2;
t168 = qJ(3) + qJ(4) + qJ(5);
t165 = sin(t168);
t160 = 0.1e1 / t165 ^ 2;
t166 = cos(t168);
t163 = t166 ^ 2;
t217 = t160 * t163;
t158 = t170 * t217 + 0.1e1;
t156 = 0.1e1 / t158;
t159 = 0.1e1 / t165;
t172 = sin(qJ(1));
t206 = qJD(1) * t172;
t198 = t166 * t206;
t167 = qJD(3) + qJD(4) + qJD(5);
t212 = t167 * t174;
t199 = t160 * t212;
t130 = ((t165 * t212 + t198) * t159 + t163 * t199) * t156;
t233 = -t130 + t212;
t209 = t174 * t166;
t155 = atan2(-t209, t165);
t149 = sin(t155);
t150 = cos(t155);
t140 = -t149 * t209 + t150 * t165;
t137 = 0.1e1 / t140;
t171 = sin(qJ(6));
t208 = t174 * t171;
t173 = cos(qJ(6));
t210 = t172 * t173;
t152 = t165 * t210 + t208;
t146 = 0.1e1 / t152;
t138 = 0.1e1 / t140 ^ 2;
t147 = 0.1e1 / t152 ^ 2;
t169 = t172 ^ 2;
t216 = t163 * t169;
t133 = t138 * t216 + 0.1e1;
t205 = qJD(1) * t174;
t190 = t163 * t172 * t205;
t214 = t166 * t167;
t221 = t150 * t166;
t229 = t130 * t174;
t125 = (t167 - t229) * t221 + (t233 * t165 + t198) * t149;
t231 = t125 * t137 * t138;
t232 = (-t216 * t231 + (-t165 * t169 * t214 + t190) * t138) / t133 ^ 2;
t193 = qJD(1) * t165 + qJD(6);
t186 = t193 * t174;
t194 = qJD(6) * t165 + qJD(1);
t187 = t194 * t173;
t213 = t167 * t172;
t135 = t172 * t187 + (t166 * t213 + t186) * t171;
t207 = t174 * t173;
t211 = t172 * t171;
t151 = t165 * t211 - t207;
t145 = t151 ^ 2;
t144 = t145 * t147 + 0.1e1;
t223 = t147 * t151;
t188 = t194 * t171;
t136 = t173 * t186 + (t173 * t214 - t188) * t172;
t227 = t136 * t146 * t147;
t230 = (t135 * t223 - t145 * t227) / t144 ^ 2;
t162 = t166 * t163;
t218 = t159 * t166;
t184 = t167 * (-t159 * t160 * t162 - t218);
t228 = (-t160 * t190 + t170 * t184) / t158 ^ 2;
t226 = t138 * t166;
t225 = t138 * t172;
t224 = t146 * t171;
t222 = t149 * t174;
t220 = t151 * t173;
t219 = t159 * t163;
t215 = t165 * t167;
t204 = -0.2e1 * t231;
t203 = 0.2e1 * t230;
t202 = t166 * t232;
t201 = t166 * t228;
t200 = t166 * t225;
t197 = 0.1e1 + t217;
t196 = t228 * t234;
t195 = 0.2e1 * t151 * t227;
t192 = t150 * t156 * t219;
t191 = (-t156 + 0.1e1) * t166 * t149;
t189 = t197 * t172;
t185 = t147 * t220 - t224;
t183 = t185 * t172;
t182 = t167 * t209 - t193 * t172;
t154 = t165 * t207 - t211;
t153 = t165 * t208 + t210;
t142 = 0.1e1 / t144;
t141 = t197 * t174 * t156;
t131 = 0.1e1 / t133;
t129 = (-t174 * t192 + t191) * t172;
t127 = t165 * t222 + t221 + (-t149 * t165 - t150 * t209) * t141;
t126 = -t197 * t196 + (-qJD(1) * t189 + t184 * t234) * t156;
t123 = t166 * t183 * t203 + (t183 * t215 + (-t185 * t205 + ((qJD(6) * t146 + t195) * t173 + (-t135 * t173 + (qJD(6) * t151 - t136) * t171) * t147) * t172) * t166) * t142;
t122 = 0.2e1 * (-t127 * t226 - t137 * t165) * t172 * t232 + ((t137 * t205 + (-t127 * t167 - t125) * t225) * t165 + (t137 * t213 + (-t126 * t150 * t174 + t233 * t149 + (t130 * t222 - t149 * t167 + t150 * t206) * t141) * t200 + (t138 * t205 + t172 * t204) * t127 + ((-t126 - t206) * t149 + ((t141 * t174 - 0.1e1) * t167 + (-t141 + t174) * t130) * t150) * t165 * t225) * t166) * t131;
t1 = [-0.2e1 * t172 * t159 * t201 + (-t167 * t189 + t205 * t218) * t156, 0, t126, t126, t126, 0; (0.2e1 * t137 * t202 + (t137 * t215 + (qJD(1) * t129 + t125) * t226) * t131) * t174 + (-0.2e1 * t138 * t202 * t129 + (((0.2e1 * t201 - t215 + (t219 * t229 + t215) * t156) * t149 + (t196 * t219 + t130 * t166 + (t162 * t199 + (-t130 + 0.2e1 * t212) * t166) * t156) * t150) * t200 + (-t138 * t215 + t166 * t204) * t129 + (t137 + ((t169 - t170) * t192 + t174 * t191) * t138) * t166 * qJD(1)) * t131) * t172, 0, t122, t122, t122, 0; (-t146 * t153 + t154 * t223) * t203 + (t154 * t195 + t174 * t146 * t187 + t182 * t224 + (t174 * t151 * t188 - t154 * t135 - t153 * t136 - t182 * t220) * t147) * t142, 0, t123, t123, t123, -0.2e1 * t230 + 0.2e1 * (t135 * t147 * t142 + (-t142 * t227 - t147 * t230) * t151) * t151;];
JaD_rot  = t1;
