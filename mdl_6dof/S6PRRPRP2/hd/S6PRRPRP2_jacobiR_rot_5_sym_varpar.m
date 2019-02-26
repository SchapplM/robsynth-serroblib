% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRP2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:47
% EndTime: 2019-02-26 20:01:48
% DurationCPUTime: 0.07s
% Computational Cost: add. (93->28), mult. (163->65), div. (0->0), fcn. (238->10), ass. (0->31)
t119 = qJ(3) + pkin(11);
t118 = cos(t119);
t124 = sin(qJ(5));
t136 = t118 * t124;
t126 = cos(qJ(5));
t135 = t118 * t126;
t120 = sin(pkin(10));
t121 = sin(pkin(6));
t134 = t120 * t121;
t122 = cos(pkin(10));
t133 = t121 * t122;
t125 = sin(qJ(2));
t132 = t121 * t125;
t123 = cos(pkin(6));
t131 = t123 * t125;
t127 = cos(qJ(2));
t130 = t123 * t127;
t129 = t124 * t127;
t128 = t126 * t127;
t117 = sin(t119);
t115 = -t120 * t131 + t122 * t127;
t114 = t120 * t130 + t122 * t125;
t113 = t120 * t127 + t122 * t131;
t112 = t120 * t125 - t122 * t130;
t111 = t123 * t117 + t118 * t132;
t110 = -t117 * t132 + t123 * t118;
t109 = t115 * t118 + t117 * t134;
t108 = -t115 * t117 + t118 * t134;
t107 = t113 * t118 - t117 * t133;
t106 = -t113 * t117 - t118 * t133;
t1 = [0, -t114 * t135 + t115 * t124, t108 * t126, 0, -t109 * t124 + t114 * t126, 0; 0, -t112 * t135 + t113 * t124, t106 * t126, 0, -t107 * t124 + t112 * t126, 0; 0 (t118 * t128 + t124 * t125) * t121, t110 * t126, 0, -t111 * t124 - t121 * t128, 0; 0, t114 * t136 + t115 * t126, -t108 * t124, 0, -t109 * t126 - t114 * t124, 0; 0, t112 * t136 + t113 * t126, -t106 * t124, 0, -t107 * t126 - t112 * t124, 0; 0 (-t118 * t129 + t125 * t126) * t121, -t110 * t124, 0, -t111 * t126 + t121 * t129, 0; 0, -t114 * t117, t109, 0, 0, 0; 0, -t112 * t117, t107, 0, 0, 0; 0, t121 * t127 * t117, t111, 0, 0, 0;];
JR_rot  = t1;
