% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:42
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:42:43
% EndTime: 2019-02-22 09:42:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (90->28), mult. (163->65), div. (0->0), fcn. (238->10), ass. (0->31)
t121 = qJ(3) + pkin(11);
t119 = sin(t121);
t126 = sin(qJ(6));
t138 = t119 * t126;
t128 = cos(qJ(6));
t137 = t119 * t128;
t122 = sin(pkin(10));
t123 = sin(pkin(6));
t136 = t122 * t123;
t124 = cos(pkin(10));
t135 = t123 * t124;
t127 = sin(qJ(2));
t134 = t123 * t127;
t125 = cos(pkin(6));
t133 = t125 * t127;
t129 = cos(qJ(2));
t132 = t125 * t129;
t131 = t126 * t129;
t130 = t128 * t129;
t120 = cos(t121);
t117 = -t122 * t133 + t124 * t129;
t116 = t122 * t132 + t124 * t127;
t115 = t122 * t129 + t124 * t133;
t114 = t122 * t127 - t124 * t132;
t113 = t125 * t119 + t120 * t134;
t112 = t119 * t134 - t125 * t120;
t111 = t117 * t120 + t119 * t136;
t110 = t117 * t119 - t120 * t136;
t109 = t115 * t120 - t119 * t135;
t108 = t115 * t119 + t120 * t135;
t1 = [0, -t116 * t138 + t117 * t128, t111 * t126, 0, 0, t110 * t128 - t116 * t126; 0, -t114 * t138 + t115 * t128, t109 * t126, 0, 0, t108 * t128 - t114 * t126; 0 (t119 * t131 + t127 * t128) * t123, t113 * t126, 0, 0, t112 * t128 + t123 * t131; 0, -t116 * t137 - t117 * t126, t111 * t128, 0, 0, -t110 * t126 - t116 * t128; 0, -t114 * t137 - t115 * t126, t109 * t128, 0, 0, -t108 * t126 - t114 * t128; 0 (t119 * t130 - t126 * t127) * t123, t113 * t128, 0, 0, -t112 * t126 + t123 * t130; 0, -t116 * t120, -t110, 0, 0, 0; 0, -t114 * t120, -t108, 0, 0, 0; 0, t123 * t129 * t120, -t112, 0, 0, 0;];
JR_rot  = t1;
