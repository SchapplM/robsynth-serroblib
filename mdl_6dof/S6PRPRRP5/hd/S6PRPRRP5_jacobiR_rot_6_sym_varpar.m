% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:36
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRP5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:36:23
% EndTime: 2019-02-22 09:36:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (57->28), mult. (163->64), div. (0->0), fcn. (238->10), ass. (0->30)
t118 = sin(pkin(6));
t122 = sin(qJ(4));
t135 = t118 * t122;
t125 = cos(qJ(4));
t134 = t118 * t125;
t126 = cos(qJ(2));
t133 = t118 * t126;
t120 = cos(pkin(6));
t123 = sin(qJ(2));
t132 = t120 * t123;
t131 = t120 * t126;
t121 = sin(qJ(5));
t130 = t121 * t122;
t129 = t121 * t123;
t124 = cos(qJ(5));
t128 = t122 * t124;
t127 = t123 * t124;
t119 = cos(pkin(10));
t117 = sin(pkin(10));
t115 = t120 * t125 - t122 * t133;
t114 = -t120 * t122 - t125 * t133;
t113 = -t117 * t132 + t119 * t126;
t112 = t117 * t131 + t119 * t123;
t111 = t117 * t126 + t119 * t132;
t110 = t117 * t123 - t119 * t131;
t109 = t110 * t122 - t119 * t134;
t108 = t110 * t125 + t119 * t135;
t107 = t112 * t122 + t117 * t134;
t106 = t112 * t125 - t117 * t135;
t1 = [0, -t112 * t121 + t113 * t128, 0, t106 * t124, -t107 * t121 + t113 * t124, 0; 0, -t110 * t121 + t111 * t128, 0, t108 * t124, -t109 * t121 + t111 * t124, 0; 0 (t121 * t126 + t122 * t127) * t118, 0, t114 * t124, -t115 * t121 + t118 * t127, 0; 0, -t112 * t124 - t113 * t130, 0, -t106 * t121, -t107 * t124 - t113 * t121, 0; 0, -t110 * t124 - t111 * t130, 0, -t108 * t121, -t109 * t124 - t111 * t121, 0; 0 (-t122 * t129 + t124 * t126) * t118, 0, -t114 * t121, -t115 * t124 - t118 * t129, 0; 0, -t113 * t125, 0, t107, 0, 0; 0, -t111 * t125, 0, t109, 0, 0; 0, -t123 * t134, 0, t115, 0, 0;];
JR_rot  = t1;
