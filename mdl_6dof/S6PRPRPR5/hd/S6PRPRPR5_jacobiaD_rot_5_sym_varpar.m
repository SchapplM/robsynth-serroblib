% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRPR5_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:45
% EndTime: 2019-02-26 19:48:46
% DurationCPUTime: 0.75s
% Computational Cost: add. (4935->90), mult. (7372->193), div. (524->12), fcn. (9540->11), ass. (0->88)
t173 = sin(pkin(10));
t175 = cos(pkin(10));
t177 = sin(qJ(2));
t176 = cos(pkin(6));
t178 = cos(qJ(2));
t199 = t176 * t178;
t162 = -t173 * t177 + t175 * t199;
t155 = t162 * qJD(2);
t200 = t176 * t177;
t163 = t173 * t178 + t175 * t200;
t172 = pkin(11) + qJ(4);
t170 = sin(t172);
t174 = sin(pkin(6));
t203 = t174 * t175;
t191 = t170 * t203;
t171 = cos(t172);
t196 = qJD(4) * t171;
t127 = -qJD(4) * t191 + t155 * t170 + t163 * t196;
t145 = t163 * t170 + t171 * t203;
t142 = t145 ^ 2;
t202 = t174 * t177;
t153 = t170 * t202 - t171 * t176;
t151 = 0.1e1 / t153 ^ 2;
t135 = t142 * t151 + 0.1e1;
t133 = 0.1e1 / t135;
t154 = t170 * t176 + t171 * t202;
t197 = qJD(2) * t178;
t190 = t174 * t197;
t140 = qJD(4) * t154 + t170 * t190;
t150 = 0.1e1 / t153;
t207 = t145 * t151;
t115 = (-t127 * t150 + t140 * t207) * t133;
t136 = atan2(-t145, t153);
t131 = sin(t136);
t132 = cos(t136);
t188 = -t131 * t153 - t132 * t145;
t112 = t115 * t188 - t127 * t131 + t132 * t140;
t126 = -t131 * t145 + t132 * t153;
t123 = 0.1e1 / t126;
t124 = 0.1e1 / t126 ^ 2;
t216 = t112 * t123 * t124;
t192 = t173 * t200;
t165 = t175 * t178 - t192;
t204 = t173 * t174;
t186 = -t165 * t170 + t171 * t204;
t215 = -0.2e1 * t186 * t216;
t201 = t174 * t178;
t185 = -t150 * t162 + t201 * t207;
t214 = t170 * t185;
t208 = t140 * t150 * t151;
t213 = -0.2e1 * (t127 * t207 - t142 * t208) / t135 ^ 2;
t164 = t173 * t199 + t175 * t177;
t159 = 0.1e1 / t164;
t160 = 0.1e1 / t164 ^ 2;
t212 = t124 * t186;
t149 = t165 * t171 + t170 * t204;
t157 = t164 * qJD(2);
t129 = qJD(4) * t149 - t157 * t170;
t211 = t129 * t124;
t210 = t131 * t186;
t209 = t132 * t186;
t206 = t149 * t165;
t205 = t164 * t170;
t198 = qJD(2) * t177;
t143 = t186 ^ 2;
t121 = t124 * t143 + 0.1e1;
t195 = 0.2e1 * (-t143 * t216 - t186 * t211) / t121 ^ 2;
t130 = qJD(4) * t186 - t157 * t171;
t144 = t149 ^ 2;
t139 = t144 * t160 + 0.1e1;
t158 = -qJD(2) * t192 + t175 * t197;
t161 = t159 * t160;
t194 = 0.2e1 * (t130 * t149 * t160 - t144 * t158 * t161) / t139 ^ 2;
t189 = -0.2e1 * t145 * t208;
t147 = t163 * t171 - t191;
t187 = -t147 * t150 + t154 * t207;
t156 = t163 * qJD(2);
t141 = -qJD(4) * t153 + t171 * t190;
t137 = 0.1e1 / t139;
t128 = -qJD(4) * t145 + t155 * t171;
t118 = 0.1e1 / t121;
t117 = t133 * t214;
t116 = t187 * t133;
t114 = (-t131 * t162 + t132 * t201) * t170 + t188 * t117;
t113 = t116 * t188 - t131 * t147 + t132 * t154;
t111 = t187 * t213 + (t154 * t189 - t128 * t150 + (t127 * t154 + t140 * t147 + t141 * t145) * t151) * t133;
t109 = t213 * t214 + (t185 * t196 + (t189 * t201 + t150 * t156 + (t140 * t162 + (t127 * t178 - t145 * t198) * t174) * t151) * t170) * t133;
t1 = [0, t109, 0, t111, 0, 0; 0 (-t114 * t212 + t123 * t205) * t195 + ((-t158 * t170 - t164 * t196) * t123 + (-t211 + t215) * t114 + (t205 * t112 + (-t109 * t145 - t117 * t127 + (-t170 * t198 + t178 * t196) * t174 + (-t117 * t153 - t162 * t170) * t115) * t209 + (-t162 * t196 - t109 * t153 - t117 * t140 + t156 * t170 + (t117 * t145 - t170 * t201) * t115) * t210) * t124) * t118, 0 (-t113 * t212 - t123 * t149) * t195 + (t113 * t215 + t130 * t123 + (-t149 * t112 - t113 * t129 + (-t111 * t145 - t116 * t127 + t141 + (-t116 * t153 - t147) * t115) * t209 + (-t111 * t153 - t116 * t140 - t128 + (t116 * t145 - t154) * t115) * t210) * t124) * t118, 0, 0; 0 (t159 * t164 * t171 + t160 * t206) * t194 + (qJD(4) * t159 * t205 + (-t130 * t165 + t149 * t157) * t160 + (0.2e1 * t161 * t206 + (t160 * t164 - t159) * t171) * t158) * t137, 0, -t186 * t159 * t194 + (-t158 * t160 * t186 - t129 * t159) * t137, 0, 0;];
JaD_rot  = t1;
