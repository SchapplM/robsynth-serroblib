% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR13_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:50
% EndTime: 2019-02-26 21:44:50
% DurationCPUTime: 0.07s
% Computational Cost: add. (57->27), mult. (163->58), div. (0->0), fcn. (238->10), ass. (0->30)
t112 = sin(pkin(11));
t116 = sin(qJ(4));
t131 = t112 * t116;
t113 = sin(pkin(6));
t119 = cos(qJ(4));
t130 = t113 * t119;
t120 = cos(qJ(2));
t129 = t113 * t120;
t121 = cos(qJ(1));
t128 = t113 * t121;
t114 = cos(pkin(11));
t127 = t114 * t116;
t117 = sin(qJ(2));
t126 = t116 * t117;
t118 = sin(qJ(1));
t125 = t118 * t117;
t124 = t118 * t120;
t123 = t121 * t117;
t122 = t121 * t120;
t115 = cos(pkin(6));
t106 = -t115 * t122 + t125;
t102 = t106 * t119 + t116 * t128;
t103 = -t106 * t116 + t119 * t128;
t109 = -t115 * t125 + t122;
t108 = t115 * t124 + t123;
t107 = t115 * t123 + t124;
t105 = -t115 * t116 - t119 * t129;
t101 = t108 * t116 + t118 * t130;
t100 = t118 * t113 * t116 - t108 * t119;
t1 = [t103 * t114 - t107 * t112, -t108 * t112 + t109 * t127, 0, -t100 * t114, 0, 0; t101 * t114 + t109 * t112, -t106 * t112 + t107 * t127, 0, t102 * t114, 0, 0; 0 (t112 * t120 + t114 * t126) * t113, 0, t105 * t114, 0, 0; -t103 * t112 - t107 * t114, -t108 * t114 - t109 * t131, 0, t100 * t112, 0, 0; -t101 * t112 + t109 * t114, -t106 * t114 - t107 * t131, 0, -t102 * t112, 0, 0; 0 (-t112 * t126 + t114 * t120) * t113, 0, -t105 * t112, 0, 0; t102, -t109 * t119, 0, t101, 0, 0; t100, -t107 * t119, 0, -t103, 0, 0; 0, -t117 * t130, 0, t115 * t119 - t116 * t129, 0, 0;];
JR_rot  = t1;
