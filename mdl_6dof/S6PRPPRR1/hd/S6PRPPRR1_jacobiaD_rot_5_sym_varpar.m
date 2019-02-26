% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:36
% DurationCPUTime: 0.52s
% Computational Cost: add. (1932->59), mult. (5333->135), div. (281->12), fcn. (6885->13), ass. (0->72)
t188 = pkin(12) + qJ(5);
t186 = sin(t188);
t187 = cos(t188);
t194 = cos(pkin(6));
t189 = sin(pkin(11));
t192 = cos(pkin(11));
t195 = sin(qJ(2));
t196 = cos(qJ(2));
t208 = t189 * t196 + t195 * t192;
t178 = t208 * t194;
t181 = t195 * t189 - t192 * t196;
t190 = sin(pkin(10));
t193 = cos(pkin(10));
t209 = -t178 * t190 - t181 * t193;
t191 = sin(pkin(6));
t212 = t190 * t191;
t205 = -t186 * t209 + t187 * t212;
t222 = t205 * qJD(5);
t179 = t181 * qJD(2);
t204 = t181 * t194;
t163 = -t190 * t208 - t193 * t204;
t176 = t181 * t191;
t149 = atan2(t163, t176);
t144 = sin(t149);
t145 = cos(t149);
t138 = t144 * t163 + t145 * t176;
t135 = 0.1e1 / t138;
t155 = t186 * t212 + t187 * t209;
t151 = 0.1e1 / t155;
t173 = 0.1e1 / t176;
t136 = 0.1e1 / t138 ^ 2;
t152 = 0.1e1 / t155 ^ 2;
t174 = 0.1e1 / t176 ^ 2;
t160 = t163 ^ 2;
t148 = t160 * t174 + 0.1e1;
t146 = 0.1e1 / t148;
t203 = qJD(2) * t178;
t156 = t190 * t179 - t193 * t203;
t177 = t208 * t191;
t170 = qJD(2) * t177;
t216 = t163 * t174;
t129 = (t156 * t173 - t170 * t216) * t146;
t210 = -t144 * t176 + t145 * t163;
t126 = t210 * t129 + t144 * t156 + t145 * t170;
t221 = t126 * t135 * t136;
t150 = t205 ^ 2;
t141 = t150 * t152 + 0.1e1;
t172 = t194 * t179;
t180 = t208 * qJD(2);
t159 = t172 * t190 - t180 * t193;
t142 = t155 * qJD(5) + t159 * t186;
t217 = t152 * t205;
t143 = t159 * t187 + t222;
t218 = t143 * t151 * t152;
t220 = (-t142 * t217 - t150 * t218) / t141 ^ 2;
t165 = t190 * t204 - t193 * t208;
t219 = t136 * t165;
t215 = t163 * t177;
t214 = t170 * t173 * t174;
t207 = -t151 * t186 - t187 * t217;
t162 = -t178 * t193 + t181 * t190;
t206 = -t162 * t173 + t174 * t215;
t171 = t191 * t179;
t161 = t165 ^ 2;
t158 = t179 * t193 + t190 * t203;
t157 = t172 * t193 + t180 * t190;
t139 = 0.1e1 / t141;
t133 = t136 * t161 + 0.1e1;
t130 = t206 * t146;
t127 = -t210 * t130 + t144 * t162 + t145 * t177;
t125 = 0.2e1 * t206 / t148 ^ 2 * (t156 * t216 - t160 * t214) + (0.2e1 * t214 * t215 + t157 * t173 + (-t156 * t177 - t162 * t170 + t163 * t171) * t174) * t146;
t1 = [0, t125, 0, 0, 0, 0; 0, 0.2e1 * (-t127 * t219 - t135 * t209) / t133 ^ 2 * (t158 * t219 - t161 * t221) + (t159 * t135 + (-t126 * t209 + t127 * t158) * t136 + (-0.2e1 * t127 * t221 + ((t125 * t163 - t130 * t156 - t171 + (t130 * t176 + t162) * t129) * t145 + (-t125 * t176 + t130 * t170 + t157 + (t130 * t163 - t177) * t129) * t144) * t136) * t165) / t133, 0, 0, 0, 0; 0, 0.2e1 * t207 * t165 * t220 + (-t207 * t158 + ((qJD(5) * t151 - 0.2e1 * t205 * t218) * t187 + (-t142 * t187 + (-t143 - t222) * t186) * t152) * t165) * t139, 0, 0, -0.2e1 * t220 - 0.2e1 * (t139 * t142 * t152 - (-t139 * t218 - t152 * t220) * t205) * t205, 0;];
JaD_rot  = t1;
