% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRPRR4_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:52
% EndTime: 2019-02-26 20:05:53
% DurationCPUTime: 0.66s
% Computational Cost: add. (1361->69), mult. (4214->164), div. (442->14), fcn. (5379->13), ass. (0->80)
t235 = -qJD(3) + qJD(5);
t192 = sin(qJ(3));
t195 = cos(qJ(3));
t187 = sin(pkin(11));
t189 = cos(pkin(11));
t196 = cos(qJ(2));
t190 = cos(pkin(6));
t193 = sin(qJ(2));
t221 = t190 * t193;
t203 = t187 * t221 - t189 * t196;
t188 = sin(pkin(6));
t223 = t187 * t188;
t162 = t192 * t223 - t195 * t203;
t191 = sin(qJ(5));
t194 = cos(qJ(5));
t205 = t192 * t203 + t195 * t223;
t155 = t162 * t194 - t191 * t205;
t220 = t190 * t196;
t204 = t187 * t220 + t189 * t193;
t171 = t204 * qJD(2);
t156 = t162 * qJD(3) - t171 * t192;
t157 = t205 * qJD(3) - t171 * t195;
t133 = t155 * qJD(5) - t156 * t194 + t157 * t191;
t150 = 0.1e1 / t155 ^ 2;
t210 = -t162 * t191 - t194 * t205;
t227 = t150 * t210;
t234 = t133 * t227;
t233 = t235 * t195;
t232 = t235 * t192;
t134 = t210 * qJD(5) + t156 * t191 + t157 * t194;
t175 = t187 * t193 - t189 * t220;
t222 = t188 * t196;
t167 = atan2(t175, t222);
t163 = sin(t167);
t164 = cos(t167);
t146 = t163 * t175 + t164 * t222;
t143 = 0.1e1 / t146;
t149 = 0.1e1 / t155;
t184 = 0.1e1 / t196;
t144 = 0.1e1 / t146 ^ 2;
t185 = 0.1e1 / t196 ^ 2;
t148 = t210 ^ 2;
t137 = t148 * t150 + 0.1e1;
t151 = t149 * t150;
t229 = t134 * t151;
t231 = (-t148 * t229 - t234) / t137 ^ 2;
t176 = t187 * t196 + t189 * t221;
t170 = t176 * qJD(2);
t224 = t185 * t193;
t212 = t175 * t224;
t173 = t175 ^ 2;
t183 = 0.1e1 / t188 ^ 2;
t168 = t173 * t183 * t185 + 0.1e1;
t165 = 0.1e1 / t168;
t182 = 0.1e1 / t188;
t225 = t165 * t182;
t138 = (qJD(2) * t212 + t170 * t184) * t225;
t207 = -t163 * t222 + t164 * t175;
t213 = t164 * t188 * t193;
t131 = -qJD(2) * t213 + t207 * t138 + t163 * t170;
t230 = t131 * t143 * t144;
t228 = t144 * t204;
t208 = t191 * t192 + t194 * t195;
t159 = t208 * t204;
t226 = t150 * t159;
t215 = 0.2e1 * t231;
t214 = -0.2e1 * t151 * t210;
t209 = t191 * t195 - t192 * t194;
t206 = t176 * t184 + t212;
t186 = t184 * t185;
t174 = t204 ^ 2;
t172 = t203 * qJD(2);
t169 = t175 * qJD(2);
t158 = t209 * t204;
t142 = t174 * t144 + 0.1e1;
t139 = t206 * t225;
t135 = 0.1e1 / t137;
t132 = t207 * t139 + t163 * t176 - t213;
t130 = (-0.2e1 * t206 / t168 ^ 2 * (qJD(2) * t173 * t186 * t193 + t170 * t175 * t185) * t183 + (t170 * t224 - t169 * t184 + (t176 * t224 + (0.2e1 * t186 * t193 ^ 2 + t184) * t175) * qJD(2)) * t165) * t182;
t1 = [0, t130, 0, 0, 0, 0; 0, 0.2e1 * (-t132 * t228 - t143 * t203) / t142 ^ 2 * (-t172 * t228 - t174 * t230) + (t171 * t143 + (-t131 * t203 - t132 * t172) * t144 - (0.2e1 * t132 * t230 + (-(-qJD(2) * t222 + t130 * t175 + t139 * t170 + (-t139 * t222 + t176) * t138) * t164 - (-t138 * t139 * t175 - t169 + (-t130 * t196 + (qJD(2) * t139 + t138) * t193) * t188) * t163) * t144) * t204) / t142, 0, 0, 0, 0; 0 (t149 * t158 + t210 * t226) * t215 + (t133 * t226 + (t158 * t150 - t159 * t214) * t134 + (t209 * t149 + t208 * t227) * t172 - ((t233 * t149 + t232 * t227) * t194 + (t232 * t149 - t233 * t227) * t191) * t204) * t135 (t149 * t155 + t210 * t227) * t215 + (0.2e1 * t234 + (t150 * t155 - t210 * t214 - t149) * t134) * t135, 0, -0.2e1 * t231 - 0.2e1 * (t133 * t150 * t135 - (-t135 * t229 - t150 * t231) * t210) * t210, 0;];
JaD_rot  = t1;
