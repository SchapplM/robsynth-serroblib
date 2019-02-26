% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP1
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

function JR_rot = S6PRRPRP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:11
% EndTime: 2019-02-26 20:01:11
% DurationCPUTime: 0.07s
% Computational Cost: add. (93->28), mult. (163->65), div. (0->0), fcn. (238->10), ass. (0->31)
t126 = qJ(3) + pkin(11);
t125 = cos(t126);
t131 = sin(qJ(5));
t143 = t125 * t131;
t133 = cos(qJ(5));
t142 = t125 * t133;
t127 = sin(pkin(10));
t128 = sin(pkin(6));
t141 = t127 * t128;
t129 = cos(pkin(10));
t140 = t128 * t129;
t132 = sin(qJ(2));
t139 = t128 * t132;
t130 = cos(pkin(6));
t138 = t130 * t132;
t134 = cos(qJ(2));
t137 = t130 * t134;
t136 = t131 * t134;
t135 = t133 * t134;
t124 = sin(t126);
t122 = -t127 * t138 + t129 * t134;
t121 = t127 * t137 + t129 * t132;
t120 = t127 * t134 + t129 * t138;
t119 = t127 * t132 - t129 * t137;
t118 = t130 * t124 + t125 * t139;
t117 = -t124 * t139 + t130 * t125;
t116 = t122 * t125 + t124 * t141;
t115 = -t122 * t124 + t125 * t141;
t114 = t120 * t125 - t124 * t140;
t113 = -t120 * t124 - t125 * t140;
t1 = [0, -t121 * t142 + t122 * t131, t115 * t133, 0, -t116 * t131 + t121 * t133, 0; 0, -t119 * t142 + t120 * t131, t113 * t133, 0, -t114 * t131 + t119 * t133, 0; 0 (t125 * t135 + t131 * t132) * t128, t117 * t133, 0, -t118 * t131 - t128 * t135, 0; 0, t121 * t143 + t122 * t133, -t115 * t131, 0, -t116 * t133 - t121 * t131, 0; 0, t119 * t143 + t120 * t133, -t113 * t131, 0, -t114 * t133 - t119 * t131, 0; 0 (-t125 * t136 + t132 * t133) * t128, -t117 * t131, 0, -t118 * t133 + t128 * t136, 0; 0, -t121 * t124, t116, 0, 0, 0; 0, -t119 * t124, t114, 0, 0, 0; 0, t128 * t134 * t124, t118, 0, 0, 0;];
JR_rot  = t1;
