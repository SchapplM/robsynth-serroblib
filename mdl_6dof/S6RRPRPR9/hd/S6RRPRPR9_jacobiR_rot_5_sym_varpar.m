% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:27
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR9_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:27:29
% EndTime: 2019-02-22 11:27:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (93->26), mult. (163->58), div. (0->0), fcn. (238->10), ass. (0->31)
t120 = pkin(11) + qJ(4);
t119 = cos(t120);
t121 = sin(pkin(12));
t138 = t119 * t121;
t123 = cos(pkin(12));
t137 = t119 * t123;
t127 = cos(qJ(2));
t136 = t119 * t127;
t122 = sin(pkin(6));
t125 = sin(qJ(2));
t135 = t122 * t125;
t126 = sin(qJ(1));
t134 = t122 * t126;
t128 = cos(qJ(1));
t133 = t122 * t128;
t132 = t126 * t125;
t131 = t126 * t127;
t130 = t128 * t125;
t129 = t128 * t127;
t124 = cos(pkin(6));
t114 = t124 * t130 + t131;
t118 = sin(t120);
t108 = -t114 * t118 - t119 * t133;
t109 = -t114 * t119 + t118 * t133;
t116 = -t124 * t132 + t129;
t115 = t124 * t131 + t130;
t113 = t124 * t129 - t132;
t112 = -t118 * t135 + t124 * t119;
t111 = t116 * t119 + t118 * t134;
t110 = t116 * t118 - t119 * t134;
t1 = [t109 * t123 + t113 * t121, -t115 * t137 + t116 * t121, 0, -t110 * t123, 0, 0; t111 * t123 + t115 * t121, t113 * t137 + t114 * t121, 0, t108 * t123, 0, 0; 0 (t121 * t125 + t123 * t136) * t122, 0, t112 * t123, 0, 0; -t109 * t121 + t113 * t123, t115 * t138 + t116 * t123, 0, t110 * t121, 0, 0; -t111 * t121 + t115 * t123, -t113 * t138 + t114 * t123, 0, -t108 * t121, 0, 0; 0 (-t121 * t136 + t123 * t125) * t122, 0, -t112 * t121, 0, 0; t108, -t115 * t118, 0, t111, 0, 0; t110, t113 * t118, 0, -t109, 0, 0; 0, t122 * t127 * t118, 0, t124 * t118 + t119 * t135, 0, 0;];
JR_rot  = t1;
