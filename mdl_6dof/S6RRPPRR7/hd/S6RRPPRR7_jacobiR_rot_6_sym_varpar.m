% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:57
% EndTime: 2019-02-26 21:31:57
% DurationCPUTime: 0.13s
% Computational Cost: add. (74->27), mult. (219->61), div. (0->0), fcn. (320->10), ass. (0->35)
t127 = cos(pkin(6));
t134 = cos(qJ(2));
t135 = cos(qJ(1));
t136 = t135 * t134;
t130 = sin(qJ(2));
t131 = sin(qJ(1));
t140 = t131 * t130;
t120 = -t127 * t136 + t140;
t129 = sin(qJ(5));
t133 = cos(qJ(5));
t126 = sin(pkin(6));
t142 = t126 * t135;
t114 = t120 * t133 + t129 * t142;
t137 = t135 * t130;
t139 = t131 * t134;
t121 = t127 * t137 + t139;
t128 = sin(qJ(6));
t132 = cos(qJ(6));
t149 = t114 * t128 - t121 * t132;
t148 = -t114 * t132 - t121 * t128;
t145 = t126 * t130;
t144 = t126 * t131;
t143 = t126 * t134;
t141 = t128 * t133;
t138 = t132 * t133;
t113 = -t120 * t129 + t133 * t142;
t123 = -t127 * t140 + t136;
t122 = t127 * t139 + t137;
t119 = -t127 * t129 - t133 * t143;
t118 = -t127 * t133 + t129 * t143;
t117 = t122 * t133 - t129 * t144;
t116 = t122 * t129 + t133 * t144;
t112 = t117 * t132 + t123 * t128;
t111 = -t117 * t128 + t123 * t132;
t1 = [t148, -t122 * t128 + t123 * t138, 0, 0, -t116 * t132, t111; t112, -t120 * t128 + t121 * t138, 0, 0, t113 * t132, -t149; 0 (t128 * t134 + t130 * t138) * t126, 0, 0, t118 * t132, -t119 * t128 + t132 * t145; t149, -t122 * t132 - t123 * t141, 0, 0, t116 * t128, -t112; t111, -t120 * t132 - t121 * t141, 0, 0, -t113 * t128, t148; 0 (-t130 * t141 + t132 * t134) * t126, 0, 0, -t118 * t128, -t119 * t132 - t128 * t145; t113, t123 * t129, 0, 0, t117, 0; t116, t121 * t129, 0, 0, t114, 0; 0, t129 * t145, 0, 0, t119, 0;];
JR_rot  = t1;
