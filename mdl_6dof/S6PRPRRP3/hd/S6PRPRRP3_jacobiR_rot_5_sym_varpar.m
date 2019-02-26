% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRP3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:34
% EndTime: 2019-02-26 19:51:34
% DurationCPUTime: 0.07s
% Computational Cost: add. (93->28), mult. (163->65), div. (0->0), fcn. (238->10), ass. (0->31)
t115 = pkin(11) + qJ(4);
t114 = cos(t115);
t120 = sin(qJ(5));
t132 = t114 * t120;
t122 = cos(qJ(5));
t131 = t114 * t122;
t116 = sin(pkin(10));
t117 = sin(pkin(6));
t130 = t116 * t117;
t118 = cos(pkin(10));
t129 = t117 * t118;
t121 = sin(qJ(2));
t128 = t117 * t121;
t119 = cos(pkin(6));
t127 = t119 * t121;
t123 = cos(qJ(2));
t126 = t119 * t123;
t125 = t120 * t123;
t124 = t122 * t123;
t113 = sin(t115);
t111 = -t116 * t127 + t118 * t123;
t110 = t116 * t126 + t118 * t121;
t109 = t116 * t123 + t118 * t127;
t108 = t116 * t121 - t118 * t126;
t107 = t119 * t113 + t114 * t128;
t106 = -t113 * t128 + t119 * t114;
t105 = t111 * t114 + t113 * t130;
t104 = -t111 * t113 + t114 * t130;
t103 = t109 * t114 - t113 * t129;
t102 = -t109 * t113 - t114 * t129;
t1 = [0, -t110 * t131 + t111 * t120, 0, t104 * t122, -t105 * t120 + t110 * t122, 0; 0, -t108 * t131 + t109 * t120, 0, t102 * t122, -t103 * t120 + t108 * t122, 0; 0 (t114 * t124 + t120 * t121) * t117, 0, t106 * t122, -t107 * t120 - t117 * t124, 0; 0, t110 * t132 + t111 * t122, 0, -t104 * t120, -t105 * t122 - t110 * t120, 0; 0, t108 * t132 + t109 * t122, 0, -t102 * t120, -t103 * t122 - t108 * t120, 0; 0 (-t114 * t125 + t121 * t122) * t117, 0, -t106 * t120, -t107 * t122 + t117 * t125, 0; 0, -t110 * t113, 0, t105, 0, 0; 0, -t108 * t113, 0, t103, 0, 0; 0, t117 * t123 * t113, 0, t107, 0, 0;];
JR_rot  = t1;
