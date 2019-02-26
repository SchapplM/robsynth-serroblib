% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRPRR9_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobia_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:24
% EndTime: 2019-02-26 20:53:24
% DurationCPUTime: 0.45s
% Computational Cost: add. (2181->53), mult. (6097->124), div. (85->9), fcn. (8377->19), ass. (0->73)
t130 = sin(pkin(7));
t128 = sin(pkin(13));
t132 = cos(pkin(13));
t138 = sin(qJ(3));
t142 = cos(qJ(3));
t148 = t142 * t128 + t138 * t132;
t118 = t148 * t130;
t134 = cos(pkin(7));
t120 = t148 * t134;
t135 = cos(pkin(6));
t133 = cos(pkin(12));
t143 = cos(qJ(1));
t151 = t143 * t133;
t129 = sin(pkin(12));
t139 = sin(qJ(1));
t154 = t139 * t129;
t122 = -t135 * t151 + t154;
t152 = t143 * t129;
t153 = t139 * t133;
t123 = t135 * t152 + t153;
t125 = t138 * t128 - t142 * t132;
t131 = sin(pkin(6));
t155 = t131 * t143;
t108 = t118 * t155 + t122 * t120 + t123 * t125;
t116 = -t122 * t130 + t134 * t155;
t137 = sin(qJ(5));
t141 = cos(qJ(5));
t165 = t108 * t137 - t116 * t141;
t98 = t108 * t141 + t116 * t137;
t113 = t135 * t118 + (t120 * t133 - t125 * t129) * t131;
t121 = -t131 * t133 * t130 + t135 * t134;
t103 = t113 * t137 - t121 * t141;
t94 = atan2(t165, t103);
t91 = sin(t94);
t92 = cos(t94);
t85 = t92 * t103 + t165 * t91;
t84 = 0.1e1 / t85 ^ 2;
t124 = -t135 * t154 + t151;
t147 = -t135 * t153 - t152;
t156 = t131 * t139;
t144 = t118 * t156 + t120 * t147 - t124 * t125;
t145 = -t130 * t147 + t134 * t156;
t99 = t137 * t144 - t141 * t145;
t164 = t84 * t99;
t100 = t137 * t145 + t141 * t144;
t117 = t125 * t130;
t119 = t125 * t134;
t110 = -t117 * t156 - t119 * t147 - t124 * t148;
t136 = sin(qJ(6));
t140 = cos(qJ(6));
t90 = t100 * t140 - t110 * t136;
t88 = 0.1e1 / t90 ^ 2;
t89 = t100 * t136 + t110 * t140;
t163 = t88 * t89;
t162 = t92 * t165;
t161 = t99 ^ 2 * t84;
t102 = 0.1e1 / t103 ^ 2;
t160 = t102 * t165;
t159 = t110 * t141;
t150 = t89 ^ 2 * t88 + 0.1e1;
t149 = -t103 * t91 + t162;
t146 = t117 * t155 + t122 * t119 - t123 * t148;
t112 = -t135 * t117 + (-t119 * t133 - t129 * t148) * t131;
t104 = t113 * t141 + t121 * t137;
t101 = 0.1e1 / t103;
t93 = 0.1e1 / (t102 * t165 ^ 2 + 0.1e1);
t87 = 0.1e1 / t90;
t86 = 0.1e1 / t150;
t83 = 0.1e1 / t85;
t82 = 0.1e1 / (0.1e1 + t161);
t81 = (-t101 * t146 - t112 * t160) * t93 * t137;
t80 = (t101 * t98 - t104 * t160) * t93;
t1 = [-t99 * t101 * t93, 0, t81, 0, t80, 0; (t165 * t83 - (-t91 + (-t101 * t162 + t91) * t93) * t161) * t82, 0 (t110 * t137 * t83 - (t149 * t81 + (t112 * t92 - t146 * t91) * t137) * t164) * t82, 0 (t100 * t83 - (t92 * t104 + t149 * t80 + t91 * t98) * t164) * t82, 0; ((t98 * t136 - t140 * t146) * t87 - (t136 * t146 + t98 * t140) * t163) * t86, 0 ((t136 * t159 - t140 * t144) * t87 - (t136 * t144 + t140 * t159) * t163) * t86, 0 (-t136 * t87 + t140 * t163) * t99 * t86, t150 * t86;];
Ja_rot  = t1;
