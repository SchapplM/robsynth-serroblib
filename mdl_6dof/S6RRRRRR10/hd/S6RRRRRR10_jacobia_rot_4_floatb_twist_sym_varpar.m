% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR10_jacobia_rot_4_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_rot_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_rot_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:17
% EndTime: 2018-11-23 11:27:18
% DurationCPUTime: 0.57s
% Computational Cost: add. (5275->81), mult. (5471->150), div. (86->9), fcn. (5747->29), ass. (0->92)
t204 = pkin(6) - qJ(2);
t195 = cos(t204);
t170 = pkin(6) + qJ(2);
t199 = cos(t170) / 0.2e1;
t157 = t195 / 0.2e1 + t199;
t177 = sin(qJ(2));
t178 = sin(qJ(1));
t182 = cos(qJ(1));
t140 = -t157 * t182 + t178 * t177;
t192 = sin(t204);
t198 = sin(t170) / 0.2e1;
t151 = t198 - t192 / 0.2e1;
t181 = cos(qJ(2));
t141 = t151 * t182 + t178 * t181;
t202 = pkin(7) + qJ(3);
t187 = sin(t202) / 0.2e1;
t203 = pkin(7) - qJ(3);
t191 = sin(t203);
t148 = t187 + t191 / 0.2e1;
t189 = cos(t202) / 0.2e1;
t194 = cos(t203);
t155 = t194 / 0.2e1 + t189;
t176 = sin(qJ(3));
t173 = sin(pkin(6));
t205 = t173 * t182;
t128 = t140 * t155 + t141 * t176 + t148 * t205;
t171 = sin(pkin(8));
t172 = sin(pkin(7));
t216 = cos(pkin(7));
t196 = t173 * t216;
t185 = t140 * t172 - t182 * t196;
t215 = cos(pkin(8));
t119 = -t128 * t171 - t185 * t215;
t200 = pkin(8) + qJ(4);
t217 = sin(t200) / 0.2e1;
t150 = t198 + t192 / 0.2e1;
t156 = t199 - t195 / 0.2e1;
t174 = cos(pkin(6));
t125 = -(t148 * t174 + t150 * t155 + t156 * t176) * t171 + (-t150 * t172 + t174 * t216) * t215;
t118 = atan2(t119, t125);
t115 = sin(t118);
t116 = cos(t118);
t109 = t115 * t119 + t116 * t125;
t108 = 0.1e1 / t109 ^ 2;
t143 = -t178 * t157 - t177 * t182;
t145 = t178 * t151 - t181 * t182;
t206 = t173 * t178;
t130 = t143 * t155 + t145 * t176 + t148 * t206;
t184 = t143 * t172 - t178 * t196;
t120 = t130 * t171 + t184 * t215;
t214 = t108 * t120;
t213 = t108 * t120 ^ 2;
t149 = t187 - t191 / 0.2e1;
t154 = t189 - t194 / 0.2e1;
t180 = cos(qJ(3));
t132 = -t143 * t149 + t145 * t180 + t154 * t206;
t201 = pkin(8) - qJ(4);
t190 = sin(t201);
t147 = t217 - t190 / 0.2e1;
t188 = cos(t200) / 0.2e1;
t193 = cos(t201);
t152 = t188 - t193 / 0.2e1;
t179 = cos(qJ(4));
t114 = t130 * t147 - t132 * t179 + t184 * t152;
t112 = 0.1e1 / t114 ^ 2;
t146 = t217 + t190 / 0.2e1;
t153 = t193 / 0.2e1 + t188;
t175 = sin(qJ(4));
t113 = -t130 * t153 - t132 * t175 + t184 * t146;
t212 = t112 * t113;
t211 = t113 ^ 2 * t112;
t210 = t116 * t119;
t124 = 0.1e1 / t125 ^ 2;
t209 = t119 * t124;
t207 = t145 * t172;
t197 = t172 * t215;
t186 = -t115 * t125 + t210;
t183 = t140 * t149 - t141 * t180 - t154 * t205;
t136 = -t149 * t150 + t154 * t174 + t156 * t180;
t135 = t143 * t180 + t145 * t149;
t134 = -t143 * t176 + t145 * t155;
t133 = -(-t150 * t176 + t155 * t156) * t171 - t156 * t197;
t123 = 0.1e1 / t125;
t122 = (t140 * t176 - t141 * t155) * t171 - t141 * t197;
t117 = 0.1e1 / (t119 ^ 2 * t124 + 0.1e1);
t111 = 0.1e1 / t114;
t110 = 0.1e1 / (0.1e1 + t211);
t107 = 0.1e1 / t109;
t106 = 0.1e1 / (0.1e1 + t213);
t105 = (t123 * t183 + t136 * t209) * t171 * t117;
t104 = (t122 * t123 - t133 * t209) * t117;
t1 = [t120 * t123 * t117, t104, t105, 0, 0, 0; (t119 * t107 + (t115 + (t123 * t210 - t115) * t117) * t213) * t106 ((-t134 * t171 - t145 * t197) * t107 + (t186 * t104 + t115 * t122 + t116 * t133) * t214) * t106 (-t132 * t171 * t107 + ((t115 * t183 - t116 * t136) * t171 + t186 * t105) * t214) * t106, 0, 0, 0; ((-t128 * t153 + t146 * t185 + t175 * t183) * t111 - (t128 * t147 + t152 * t185 + t179 * t183) * t212) * t110 ((-t134 * t153 + t135 * t175 + t146 * t207) * t111 - (t134 * t147 + t135 * t179 + t152 * t207) * t212) * t110 ((t130 * t175 - t132 * t153) * t111 - (t130 * t179 + t132 * t147) * t212) * t110 (t114 * t111 + t211) * t110, 0, 0;];
Ja_rot  = t1;
