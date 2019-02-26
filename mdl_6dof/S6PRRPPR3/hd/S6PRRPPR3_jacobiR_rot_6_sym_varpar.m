% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPPR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:59:24
% EndTime: 2019-02-26 19:59:24
% DurationCPUTime: 0.07s
% Computational Cost: add. (54->28), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
t109 = sin(pkin(6));
t113 = sin(qJ(3));
t125 = t109 * t113;
t116 = cos(qJ(3));
t124 = t109 * t116;
t111 = cos(pkin(6));
t114 = sin(qJ(2));
t123 = t111 * t114;
t117 = cos(qJ(2));
t122 = t111 * t117;
t112 = sin(qJ(6));
t121 = t112 * t113;
t120 = t112 * t117;
t115 = cos(qJ(6));
t119 = t113 * t115;
t118 = t115 * t117;
t110 = cos(pkin(10));
t108 = sin(pkin(10));
t106 = t111 * t113 + t114 * t124;
t105 = -t111 * t116 + t114 * t125;
t104 = -t108 * t123 + t110 * t117;
t103 = -t108 * t122 - t110 * t114;
t102 = t108 * t117 + t110 * t123;
t101 = -t108 * t114 + t110 * t122;
t100 = t104 * t116 + t108 * t125;
t99 = t104 * t113 - t108 * t124;
t98 = t102 * t116 - t110 * t125;
t97 = t102 * t113 + t110 * t124;
t1 = [0, t103 * t119 - t104 * t112, t100 * t115, 0, 0, t103 * t115 - t99 * t112; 0, t101 * t119 - t102 * t112, t98 * t115, 0, 0, t101 * t115 - t97 * t112; 0 (-t112 * t114 + t113 * t118) * t109, t106 * t115, 0, 0, -t105 * t112 + t109 * t118; 0, -t103 * t121 - t104 * t115, -t100 * t112, 0, 0, -t103 * t112 - t99 * t115; 0, -t101 * t121 - t102 * t115, -t98 * t112, 0, 0, -t101 * t112 - t97 * t115; 0 (-t113 * t120 - t114 * t115) * t109, -t106 * t112, 0, 0, -t105 * t115 - t109 * t120; 0, t103 * t116, -t99, 0, 0, 0; 0, t101 * t116, -t97, 0, 0, 0; 0, t117 * t124, -t105, 0, 0, 0;];
JR_rot  = t1;
