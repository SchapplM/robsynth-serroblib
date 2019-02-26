% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:46:27
% EndTime: 2019-02-26 19:46:28
% DurationCPUTime: 0.50s
% Computational Cost: add. (1932->59), mult. (5333->135), div. (281->12), fcn. (6885->13), ass. (0->72)
t190 = qJ(4) + pkin(12);
t188 = sin(t190);
t189 = cos(t190);
t196 = cos(pkin(6));
t191 = sin(pkin(11));
t194 = cos(pkin(11));
t197 = sin(qJ(2));
t198 = cos(qJ(2));
t210 = t191 * t198 + t197 * t194;
t180 = t210 * t196;
t183 = t197 * t191 - t194 * t198;
t192 = sin(pkin(10));
t195 = cos(pkin(10));
t211 = -t180 * t192 - t183 * t195;
t193 = sin(pkin(6));
t214 = t192 * t193;
t207 = -t188 * t211 + t189 * t214;
t224 = t207 * qJD(4);
t181 = t183 * qJD(2);
t206 = t183 * t196;
t165 = -t192 * t210 - t195 * t206;
t178 = t183 * t193;
t151 = atan2(t165, t178);
t146 = sin(t151);
t147 = cos(t151);
t140 = t146 * t165 + t147 * t178;
t137 = 0.1e1 / t140;
t157 = t188 * t214 + t189 * t211;
t153 = 0.1e1 / t157;
t175 = 0.1e1 / t178;
t138 = 0.1e1 / t140 ^ 2;
t154 = 0.1e1 / t157 ^ 2;
t176 = 0.1e1 / t178 ^ 2;
t162 = t165 ^ 2;
t150 = t162 * t176 + 0.1e1;
t148 = 0.1e1 / t150;
t205 = qJD(2) * t180;
t158 = t192 * t181 - t195 * t205;
t179 = t210 * t193;
t172 = qJD(2) * t179;
t218 = t165 * t176;
t131 = (t158 * t175 - t172 * t218) * t148;
t212 = -t146 * t178 + t147 * t165;
t128 = t212 * t131 + t146 * t158 + t147 * t172;
t223 = t128 * t137 * t138;
t152 = t207 ^ 2;
t143 = t152 * t154 + 0.1e1;
t174 = t196 * t181;
t182 = t210 * qJD(2);
t161 = t174 * t192 - t182 * t195;
t144 = t157 * qJD(4) + t161 * t188;
t219 = t154 * t207;
t145 = t161 * t189 + t224;
t220 = t145 * t153 * t154;
t222 = (-t144 * t219 - t152 * t220) / t143 ^ 2;
t167 = t192 * t206 - t195 * t210;
t221 = t138 * t167;
t217 = t165 * t179;
t216 = t172 * t175 * t176;
t209 = -t153 * t188 - t189 * t219;
t164 = -t180 * t195 + t183 * t192;
t208 = -t164 * t175 + t176 * t217;
t173 = t193 * t181;
t163 = t167 ^ 2;
t160 = t181 * t195 + t192 * t205;
t159 = t174 * t195 + t182 * t192;
t141 = 0.1e1 / t143;
t135 = t138 * t163 + 0.1e1;
t132 = t208 * t148;
t129 = -t212 * t132 + t146 * t164 + t147 * t179;
t127 = 0.2e1 * t208 / t150 ^ 2 * (t158 * t218 - t162 * t216) + (0.2e1 * t216 * t217 + t159 * t175 + (-t158 * t179 - t164 * t172 + t165 * t173) * t176) * t148;
t1 = [0, t127, 0, 0, 0, 0; 0, 0.2e1 * (-t129 * t221 - t137 * t211) / t135 ^ 2 * (t160 * t221 - t163 * t223) + (t161 * t137 + (-t128 * t211 + t129 * t160) * t138 + (-0.2e1 * t129 * t223 + ((t127 * t165 - t132 * t158 - t173 + (t132 * t178 + t164) * t131) * t147 + (-t127 * t178 + t132 * t172 + t159 + (t132 * t165 - t179) * t131) * t146) * t138) * t167) / t135, 0, 0, 0, 0; 0, 0.2e1 * t209 * t167 * t222 + (-t209 * t160 + ((qJD(4) * t153 - 0.2e1 * t207 * t220) * t189 + (-t144 * t189 + (-t145 - t224) * t188) * t154) * t167) * t141, 0, -0.2e1 * t222 - 0.2e1 * (t141 * t144 * t154 - (-t141 * t220 - t154 * t222) * t207) * t207, 0, 0;];
JaD_rot  = t1;
