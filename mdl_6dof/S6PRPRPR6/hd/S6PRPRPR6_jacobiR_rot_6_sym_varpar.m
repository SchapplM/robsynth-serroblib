% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:12
% EndTime: 2019-02-26 19:49:13
% DurationCPUTime: 0.10s
% Computational Cost: add. (87->29), mult. (163->64), div. (0->0), fcn. (238->10), ass. (0->31)
t114 = pkin(11) + qJ(6);
t112 = sin(t114);
t119 = sin(qJ(4));
t131 = t112 * t119;
t113 = cos(t114);
t130 = t113 * t119;
t116 = sin(pkin(6));
t129 = t116 * t119;
t120 = sin(qJ(2));
t128 = t116 * t120;
t121 = cos(qJ(4));
t127 = t116 * t121;
t122 = cos(qJ(2));
t126 = t116 * t122;
t118 = cos(pkin(6));
t125 = t118 * t120;
t124 = t118 * t122;
t123 = t119 * t120;
t117 = cos(pkin(10));
t115 = sin(pkin(10));
t110 = t118 * t121 - t119 * t126;
t109 = -t118 * t119 - t121 * t126;
t108 = -t115 * t125 + t117 * t122;
t107 = t115 * t124 + t117 * t120;
t106 = t115 * t122 + t117 * t125;
t105 = t115 * t120 - t117 * t124;
t104 = t105 * t119 - t117 * t127;
t103 = t105 * t121 + t117 * t129;
t102 = t107 * t119 + t115 * t127;
t101 = t107 * t121 - t115 * t129;
t1 = [0, -t107 * t112 + t108 * t130, 0, t101 * t113, 0, -t102 * t112 + t108 * t113; 0, -t105 * t112 + t106 * t130, 0, t103 * t113, 0, -t104 * t112 + t106 * t113; 0 (t112 * t122 + t113 * t123) * t116, 0, t109 * t113, 0, -t110 * t112 + t113 * t128; 0, -t107 * t113 - t108 * t131, 0, -t101 * t112, 0, -t102 * t113 - t108 * t112; 0, -t105 * t113 - t106 * t131, 0, -t103 * t112, 0, -t104 * t113 - t106 * t112; 0 (-t112 * t123 + t113 * t122) * t116, 0, -t109 * t112, 0, -t110 * t113 - t112 * t128; 0, -t108 * t121, 0, t102, 0, 0; 0, -t106 * t121, 0, t104, 0, 0; 0, -t120 * t127, 0, t110, 0, 0;];
JR_rot  = t1;
