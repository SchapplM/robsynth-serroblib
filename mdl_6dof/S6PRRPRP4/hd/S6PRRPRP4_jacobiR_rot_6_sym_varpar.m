% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRP4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:04
% EndTime: 2019-02-26 20:03:04
% DurationCPUTime: 0.06s
% Computational Cost: add. (51->27), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
t119 = sin(pkin(6));
t123 = sin(qJ(3));
t135 = t119 * t123;
t126 = cos(qJ(3));
t134 = t119 * t126;
t121 = cos(pkin(6));
t124 = sin(qJ(2));
t133 = t121 * t124;
t127 = cos(qJ(2));
t132 = t121 * t127;
t122 = sin(qJ(5));
t131 = t122 * t123;
t130 = t122 * t127;
t125 = cos(qJ(5));
t129 = t123 * t125;
t128 = t125 * t127;
t120 = cos(pkin(10));
t118 = sin(pkin(10));
t116 = t121 * t123 + t124 * t134;
t115 = -t121 * t126 + t124 * t135;
t114 = -t118 * t133 + t120 * t127;
t113 = t118 * t132 + t120 * t124;
t112 = t118 * t127 + t120 * t133;
t111 = t118 * t124 - t120 * t132;
t110 = t114 * t126 + t118 * t135;
t109 = t114 * t123 - t118 * t134;
t108 = t112 * t126 - t120 * t135;
t107 = t112 * t123 + t120 * t134;
t1 = [0, -t113 * t131 + t114 * t125, t110 * t122, 0, t109 * t125 - t113 * t122, 0; 0, -t111 * t131 + t112 * t125, t108 * t122, 0, t107 * t125 - t111 * t122, 0; 0 (t123 * t130 + t124 * t125) * t119, t116 * t122, 0, t115 * t125 + t119 * t130, 0; 0, -t113 * t129 - t114 * t122, t110 * t125, 0, -t109 * t122 - t113 * t125, 0; 0, -t111 * t129 - t112 * t122, t108 * t125, 0, -t107 * t122 - t111 * t125, 0; 0 (-t122 * t124 + t123 * t128) * t119, t116 * t125, 0, -t115 * t122 + t119 * t128, 0; 0, -t113 * t126, -t109, 0, 0, 0; 0, -t111 * t126, -t107, 0, 0, 0; 0, t127 * t134, -t115, 0, 0, 0;];
JR_rot  = t1;
