% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:42
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:42:15
% EndTime: 2019-02-22 09:42:15
% DurationCPUTime: 0.07s
% Computational Cost: add. (123->29), mult. (163->64), div. (0->0), fcn. (238->10), ass. (0->32)
t126 = pkin(12) + qJ(6);
t122 = sin(t126);
t127 = qJ(3) + pkin(11);
t125 = cos(t127);
t142 = t122 * t125;
t124 = cos(t126);
t141 = t124 * t125;
t133 = cos(qJ(2));
t140 = t125 * t133;
t128 = sin(pkin(10));
t129 = sin(pkin(6));
t139 = t128 * t129;
t130 = cos(pkin(10));
t138 = t129 * t130;
t132 = sin(qJ(2));
t137 = t129 * t132;
t136 = t129 * t133;
t131 = cos(pkin(6));
t135 = t131 * t132;
t134 = t131 * t133;
t123 = sin(t127);
t120 = -t128 * t135 + t130 * t133;
t119 = t128 * t134 + t130 * t132;
t118 = t128 * t133 + t130 * t135;
t117 = t128 * t132 - t130 * t134;
t116 = t131 * t123 + t125 * t137;
t115 = -t123 * t137 + t131 * t125;
t114 = t120 * t125 + t123 * t139;
t113 = -t120 * t123 + t125 * t139;
t112 = t118 * t125 - t123 * t138;
t111 = -t118 * t123 - t125 * t138;
t1 = [0, -t119 * t141 + t120 * t122, t113 * t124, 0, 0, -t114 * t122 + t119 * t124; 0, -t117 * t141 + t118 * t122, t111 * t124, 0, 0, -t112 * t122 + t117 * t124; 0 (t122 * t132 + t124 * t140) * t129, t115 * t124, 0, 0, -t116 * t122 - t124 * t136; 0, t119 * t142 + t120 * t124, -t113 * t122, 0, 0, -t114 * t124 - t119 * t122; 0, t117 * t142 + t118 * t124, -t111 * t122, 0, 0, -t112 * t124 - t117 * t122; 0 (-t122 * t140 + t124 * t132) * t129, -t115 * t122, 0, 0, -t116 * t124 + t122 * t136; 0, -t119 * t123, t114, 0, 0, 0; 0, -t117 * t123, t112, 0, 0, 0; 0, t123 * t136, t116, 0, 0, 0;];
JR_rot  = t1;
