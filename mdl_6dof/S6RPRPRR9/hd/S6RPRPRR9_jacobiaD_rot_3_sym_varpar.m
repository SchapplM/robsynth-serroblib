% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RPRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRR9_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiaD_rot_3_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:24
% EndTime: 2019-02-26 20:53:24
% DurationCPUTime: 0.43s
% Computational Cost: add. (994->58), mult. (3107->139), div. (132->12), fcn. (4021->13), ass. (0->80)
t183 = sin(pkin(7));
t184 = sin(pkin(6));
t185 = cos(pkin(12));
t186 = cos(pkin(7));
t187 = cos(pkin(6));
t173 = -t184 * t185 * t183 + t187 * t186;
t191 = cos(qJ(1));
t211 = t191 * t185;
t182 = sin(pkin(12));
t189 = sin(qJ(1));
t214 = t189 * t182;
t174 = -t187 * t211 + t214;
t216 = t184 * t191;
t204 = -t174 * t183 + t186 * t216;
t157 = atan2(t204, t173);
t152 = sin(t157);
t153 = cos(t157);
t138 = t152 * t204 + t153 * t173;
t135 = 0.1e1 / t138;
t188 = sin(qJ(3));
t190 = cos(qJ(3));
t200 = t187 * t214 - t211;
t212 = t191 * t182;
t213 = t189 * t185;
t201 = t187 * t213 + t212;
t217 = t184 * t189;
t209 = t183 * t217;
t202 = -t186 * t201 + t209;
t151 = t202 * t188 - t190 * t200;
t145 = 0.1e1 / t151;
t230 = t204 ^ 2;
t170 = 0.1e1 / t173;
t136 = 0.1e1 / t138 ^ 2;
t146 = 0.1e1 / t151 ^ 2;
t171 = 0.1e1 / t173 ^ 2;
t229 = -0.2e1 * t170 * t171;
t175 = -t187 * t212 - t213;
t167 = t175 * qJD(1);
t166 = t174 * qJD(1);
t210 = qJD(1) * t184;
t206 = t191 * t210;
t199 = t166 * t186 + t183 * t206;
t139 = t151 * qJD(3) + t167 * t188 - t199 * t190;
t215 = t186 * t190;
t218 = t200 * t188;
t150 = -t190 * t209 + t201 * t215 - t218;
t144 = t150 ^ 2;
t143 = t144 * t146 + 0.1e1;
t225 = t146 * t150;
t140 = t167 * t190 + t199 * t188 + (t202 * t190 + t218) * qJD(3);
t226 = t140 * t145 * t146;
t228 = (t139 * t225 - t144 * t226) / t143 ^ 2;
t168 = t201 * qJD(1);
t207 = t189 * t210;
t159 = -t168 * t183 - t186 * t207;
t156 = t230 * t171 + 0.1e1;
t154 = 0.1e1 / t156;
t198 = t152 + (t153 * t170 * t204 - t152) * t154;
t130 = t198 * t159;
t227 = t130 * t135 * t136;
t208 = t183 * t216;
t203 = t174 * t186 + t208;
t149 = t175 * t190 + t203 * t188;
t224 = t149 * t150;
t223 = t154 * t170;
t155 = 0.1e1 / t156 ^ 2;
t222 = t155 * t204;
t158 = t166 * t183 - t186 * t206;
t221 = t158 * t136;
t163 = -t183 * t201 - t186 * t217;
t220 = t159 * t163;
t219 = t175 * t188;
t205 = t183 * t207;
t169 = t200 * qJD(1);
t160 = t163 ^ 2;
t148 = -t203 * t190 + t219;
t141 = 0.1e1 / t143;
t134 = t160 * t136 + 0.1e1;
t131 = t198 * t163;
t1 = [t220 * t222 * t229 + t158 * t223, 0, 0, 0, 0, 0; 0.2e1 * (-t131 * t136 * t163 - t135 * t204) / t134 ^ 2 * (-t160 * t227 + t163 * t221) + (t159 * t135 + (-t130 * t204 + t131 * t158) * t136 + (-0.2e1 * t131 * t227 + t198 * t221 + (t152 * t171 * t222 + (0.2e1 * t223 + (t230 * t229 - t170) * t155) * t153) * t136 * t220) * t163) / t134, 0, 0, 0, 0, 0; 0.2e1 * (-t145 * t148 + t146 * t224) * t228 + ((-t168 * t215 + t169 * t188 + t190 * t205) * t145 + 0.2e1 * t224 * t226 + (-t148 * t140 - (t169 * t190 + (t168 * t186 - t205) * t188) * t150 - t149 * t139) * t146 + (t149 * t145 - (t174 * t215 + t190 * t208 - t219) * t225) * qJD(3)) * t141, 0, -0.2e1 * t228 + 0.2e1 * (t139 * t146 * t141 + (-t141 * t226 - t146 * t228) * t150) * t150, 0, 0, 0;];
JaD_rot  = t1;
