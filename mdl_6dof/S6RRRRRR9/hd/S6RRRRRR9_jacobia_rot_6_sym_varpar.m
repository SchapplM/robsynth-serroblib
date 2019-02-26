% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:05
% EndTime: 2019-02-26 22:52:07
% DurationCPUTime: 0.65s
% Computational Cost: add. (2101->66), mult. (5738->154), div. (120->9), fcn. (7872->17), ass. (0->82)
t156 = cos(pkin(6));
t162 = cos(qJ(2));
t193 = sin(qJ(1));
t172 = t193 * t162;
t159 = sin(qJ(2));
t163 = cos(qJ(1));
t176 = t163 * t159;
t146 = t156 * t176 + t172;
t158 = sin(qJ(3));
t161 = cos(qJ(3));
t173 = t193 * t159;
t175 = t163 * t162;
t145 = -t156 * t175 + t173;
t153 = sin(pkin(7));
t155 = cos(pkin(7));
t154 = sin(pkin(6));
t180 = t154 * t163;
t168 = t145 * t155 + t153 * t180;
t132 = -t146 * t161 + t168 * t158;
t142 = -t145 * t153 + t155 * t180;
t157 = sin(qJ(4));
t160 = cos(qJ(4));
t194 = t132 * t157 - t142 * t160;
t120 = t132 * t160 + t142 * t157;
t179 = t155 * t158;
t182 = t153 * t156;
t141 = t158 * t182 + (t159 * t161 + t162 * t179) * t154;
t144 = -t154 * t162 * t153 + t156 * t155;
t127 = t141 * t157 - t144 * t160;
t116 = atan2(t194, t127);
t113 = sin(t116);
t114 = cos(t116);
t107 = t113 * t194 + t114 * t127;
t106 = 0.1e1 / t107 ^ 2;
t167 = t156 * t172 + t176;
t174 = t154 * t193;
t170 = t153 * t174;
t147 = -t156 * t173 + t175;
t183 = t147 * t161;
t134 = t183 + (-t167 * t155 + t170) * t158;
t164 = t167 * t153 + t155 * t174;
t121 = t134 * t157 - t164 * t160;
t192 = t106 * t121;
t122 = t134 * t160 + t164 * t157;
t166 = t167 * t161;
t133 = t147 * t158 + t155 * t166 - t161 * t170;
t152 = qJ(5) + qJ(6);
t150 = sin(t152);
t151 = cos(t152);
t112 = t122 * t151 + t133 * t150;
t110 = 0.1e1 / t112 ^ 2;
t111 = t122 * t150 - t133 * t151;
t191 = t110 * t111;
t190 = t114 * t194;
t126 = 0.1e1 / t127 ^ 2;
t189 = t194 * t126;
t188 = t121 ^ 2 * t106;
t187 = t133 * t160;
t181 = t153 * t160;
t178 = t158 * t159;
t177 = t161 * t162;
t171 = t111 ^ 2 * t110 + 0.1e1;
t169 = -t113 * t127 + t190;
t165 = -t146 * t158 - t168 * t161;
t140 = t161 * t182 + (t155 * t177 - t178) * t154;
t137 = ((-t155 * t178 + t177) * t157 - t159 * t181) * t154;
t136 = -t147 * t179 - t166;
t135 = t155 * t183 - t167 * t158;
t128 = t141 * t160 + t144 * t157;
t125 = 0.1e1 / t127;
t124 = t147 * t153 * t157 + t136 * t160;
t123 = (-t145 * t161 - t146 * t179) * t157 - t146 * t181;
t115 = 0.1e1 / (t126 * t194 ^ 2 + 0.1e1);
t109 = 0.1e1 / t112;
t108 = 0.1e1 / t171;
t105 = 0.1e1 / t107;
t104 = 0.1e1 / (0.1e1 + t188);
t103 = (-t125 * t165 - t140 * t189) * t157 * t115;
t102 = (-t123 * t125 - t137 * t189) * t115;
t101 = (t120 * t125 - t128 * t189) * t115;
t100 = t171 * t108;
t1 = [-t121 * t125 * t115, t102, t103, t101, 0, 0; (t194 * t105 - (-t113 + (-t125 * t190 + t113) * t115) * t188) * t104 ((t136 * t157 - t147 * t181) * t105 - (t169 * t102 - t113 * t123 + t114 * t137) * t192) * t104 (-t133 * t157 * t105 - ((-t113 * t165 + t114 * t140) * t157 + t169 * t103) * t192) * t104 (t122 * t105 - (t169 * t101 + t113 * t120 + t114 * t128) * t192) * t104, 0, 0; ((t120 * t150 - t151 * t165) * t109 - (t120 * t151 + t150 * t165) * t191) * t108 ((t124 * t150 - t135 * t151) * t109 - (t124 * t151 + t135 * t150) * t191) * t108 ((-t134 * t151 - t150 * t187) * t109 - (t134 * t150 - t151 * t187) * t191) * t108 (-t150 * t109 + t151 * t191) * t121 * t108, t100, t100;];
Ja_rot  = t1;
