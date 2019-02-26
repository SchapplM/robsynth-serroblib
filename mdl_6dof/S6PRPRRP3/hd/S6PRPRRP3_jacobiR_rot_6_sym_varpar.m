% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRP3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:34
% EndTime: 2019-02-26 19:51:34
% DurationCPUTime: 0.07s
% Computational Cost: add. (93->28), mult. (163->65), div. (0->0), fcn. (238->10), ass. (0->31)
t123 = pkin(11) + qJ(4);
t122 = cos(t123);
t128 = sin(qJ(5));
t140 = t122 * t128;
t130 = cos(qJ(5));
t139 = t122 * t130;
t124 = sin(pkin(10));
t125 = sin(pkin(6));
t138 = t124 * t125;
t126 = cos(pkin(10));
t137 = t125 * t126;
t129 = sin(qJ(2));
t136 = t125 * t129;
t127 = cos(pkin(6));
t135 = t127 * t129;
t131 = cos(qJ(2));
t134 = t127 * t131;
t133 = t128 * t131;
t132 = t130 * t131;
t121 = sin(t123);
t119 = -t124 * t135 + t126 * t131;
t118 = t124 * t134 + t126 * t129;
t117 = t124 * t131 + t126 * t135;
t116 = t124 * t129 - t126 * t134;
t115 = t127 * t121 + t122 * t136;
t114 = -t121 * t136 + t127 * t122;
t113 = t119 * t122 + t121 * t138;
t112 = -t119 * t121 + t122 * t138;
t111 = t117 * t122 - t121 * t137;
t110 = -t117 * t121 - t122 * t137;
t1 = [0, -t118 * t139 + t119 * t128, 0, t112 * t130, -t113 * t128 + t118 * t130, 0; 0, -t116 * t139 + t117 * t128, 0, t110 * t130, -t111 * t128 + t116 * t130, 0; 0 (t122 * t132 + t128 * t129) * t125, 0, t114 * t130, -t115 * t128 - t125 * t132, 0; 0, t118 * t140 + t119 * t130, 0, -t112 * t128, -t113 * t130 - t118 * t128, 0; 0, t116 * t140 + t117 * t130, 0, -t110 * t128, -t111 * t130 - t116 * t128, 0; 0 (-t122 * t133 + t129 * t130) * t125, 0, -t114 * t128, -t115 * t130 + t125 * t133, 0; 0, -t118 * t121, 0, t113, 0, 0; 0, -t116 * t121, 0, t111, 0, 0; 0, t125 * t131 * t121, 0, t115, 0, 0;];
JR_rot  = t1;
