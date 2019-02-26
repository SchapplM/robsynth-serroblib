% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRP2_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:54
% EndTime: 2019-02-26 19:50:55
% DurationCPUTime: 0.49s
% Computational Cost: add. (1747->58), mult. (5333->135), div. (281->12), fcn. (6885->13), ass. (0->71)
t185 = sin(qJ(4));
t187 = cos(qJ(4));
t184 = cos(pkin(6));
t179 = sin(pkin(11));
t182 = cos(pkin(11));
t186 = sin(qJ(2));
t188 = cos(qJ(2));
t200 = t188 * t179 + t186 * t182;
t171 = t200 * t184;
t174 = t186 * t179 - t188 * t182;
t180 = sin(pkin(10));
t183 = cos(pkin(10));
t201 = -t180 * t171 - t183 * t174;
t181 = sin(pkin(6));
t205 = t180 * t181;
t197 = -t185 * t201 + t187 * t205;
t214 = qJD(4) * t197;
t172 = t174 * qJD(2);
t196 = t174 * t184;
t156 = -t180 * t200 - t183 * t196;
t169 = t174 * t181;
t142 = atan2(t156, t169);
t137 = sin(t142);
t138 = cos(t142);
t131 = t137 * t156 + t138 * t169;
t128 = 0.1e1 / t131;
t148 = t185 * t205 + t187 * t201;
t144 = 0.1e1 / t148;
t166 = 0.1e1 / t169;
t129 = 0.1e1 / t131 ^ 2;
t145 = 0.1e1 / t148 ^ 2;
t167 = 0.1e1 / t169 ^ 2;
t153 = t156 ^ 2;
t141 = t153 * t167 + 0.1e1;
t139 = 0.1e1 / t141;
t195 = qJD(2) * t171;
t149 = t180 * t172 - t183 * t195;
t170 = t200 * t181;
t163 = qJD(2) * t170;
t208 = t156 * t167;
t122 = (t149 * t166 - t163 * t208) * t139;
t202 = -t137 * t169 + t138 * t156;
t119 = t122 * t202 + t137 * t149 + t138 * t163;
t213 = t119 * t128 * t129;
t143 = t197 ^ 2;
t134 = t143 * t145 + 0.1e1;
t165 = t184 * t172;
t173 = t200 * qJD(2);
t152 = t180 * t165 - t183 * t173;
t135 = qJD(4) * t148 + t152 * t185;
t209 = t145 * t197;
t136 = t152 * t187 + t214;
t210 = t136 * t144 * t145;
t212 = (-t135 * t209 - t143 * t210) / t134 ^ 2;
t158 = t180 * t196 - t183 * t200;
t211 = t129 * t158;
t207 = t156 * t170;
t206 = t163 * t166 * t167;
t199 = -t144 * t185 - t187 * t209;
t155 = -t183 * t171 + t180 * t174;
t198 = -t155 * t166 + t167 * t207;
t164 = t181 * t172;
t154 = t158 ^ 2;
t151 = t183 * t172 + t180 * t195;
t150 = t183 * t165 + t180 * t173;
t132 = 0.1e1 / t134;
t126 = t154 * t129 + 0.1e1;
t123 = t198 * t139;
t120 = -t123 * t202 + t137 * t155 + t138 * t170;
t118 = 0.2e1 * t198 / t141 ^ 2 * (t149 * t208 - t153 * t206) + (0.2e1 * t206 * t207 + t150 * t166 + (-t149 * t170 - t155 * t163 + t156 * t164) * t167) * t139;
t1 = [0, t118, 0, 0, 0, 0; 0, 0.2e1 * (-t120 * t211 - t128 * t201) / t126 ^ 2 * (t151 * t211 - t154 * t213) + (t152 * t128 + (-t119 * t201 + t120 * t151) * t129 + (-0.2e1 * t120 * t213 + ((t118 * t156 - t123 * t149 - t164 + (t123 * t169 + t155) * t122) * t138 + (-t118 * t169 + t123 * t163 + t150 + (t123 * t156 - t170) * t122) * t137) * t129) * t158) / t126, 0, 0, 0, 0; 0, 0.2e1 * t199 * t158 * t212 + (-t199 * t151 + ((qJD(4) * t144 - 0.2e1 * t197 * t210) * t187 + (-t135 * t187 + (-t136 - t214) * t185) * t145) * t158) * t132, 0, -0.2e1 * t212 - 0.2e1 * (t132 * t135 * t145 - (-t132 * t210 - t145 * t212) * t197) * t197, 0, 0;];
JaD_rot  = t1;
