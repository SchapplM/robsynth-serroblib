% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:19
% EndTime: 2019-02-26 20:12:19
% DurationCPUTime: 0.07s
% Computational Cost: add. (84->28), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->30)
t116 = qJ(4) + pkin(12);
t114 = sin(t116);
t123 = cos(qJ(3));
t132 = t114 * t123;
t115 = cos(t116);
t131 = t115 * t123;
t118 = sin(pkin(6));
t121 = sin(qJ(3));
t130 = t118 * t121;
t129 = t118 * t123;
t124 = cos(qJ(2));
t128 = t118 * t124;
t120 = cos(pkin(6));
t122 = sin(qJ(2));
t127 = t120 * t122;
t126 = t120 * t124;
t125 = t123 * t124;
t119 = cos(pkin(11));
t117 = sin(pkin(11));
t112 = t120 * t121 + t122 * t129;
t111 = t120 * t123 - t122 * t130;
t110 = -t117 * t127 + t119 * t124;
t109 = t117 * t126 + t119 * t122;
t108 = t117 * t124 + t119 * t127;
t107 = t117 * t122 - t119 * t126;
t106 = t110 * t123 + t117 * t130;
t105 = -t110 * t121 + t117 * t129;
t104 = t108 * t123 - t119 * t130;
t103 = -t108 * t121 - t119 * t129;
t1 = [0, -t109 * t131 + t110 * t114, t105 * t115, -t106 * t114 + t109 * t115, 0, 0; 0, -t107 * t131 + t108 * t114, t103 * t115, -t104 * t114 + t107 * t115, 0, 0; 0 (t114 * t122 + t115 * t125) * t118, t111 * t115, -t112 * t114 - t115 * t128, 0, 0; 0, t109 * t132 + t110 * t115, -t105 * t114, -t106 * t115 - t109 * t114, 0, 0; 0, t107 * t132 + t108 * t115, -t103 * t114, -t104 * t115 - t107 * t114, 0, 0; 0 (-t114 * t125 + t115 * t122) * t118, -t111 * t114, -t112 * t115 + t114 * t128, 0, 0; 0, -t109 * t121, t106, 0, 0, 0; 0, -t107 * t121, t104, 0, 0, 0; 0, t121 * t128, t112, 0, 0, 0;];
JR_rot  = t1;
