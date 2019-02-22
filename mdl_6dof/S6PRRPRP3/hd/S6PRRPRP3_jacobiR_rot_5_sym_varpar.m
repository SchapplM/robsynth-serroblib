% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRP3
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
% Datum: 2019-02-22 09:46
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRP3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:46:13
% EndTime: 2019-02-22 09:46:14
% DurationCPUTime: 0.07s
% Computational Cost: add. (84->28), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->30)
t115 = pkin(11) + qJ(5);
t113 = sin(t115);
t122 = cos(qJ(3));
t131 = t113 * t122;
t114 = cos(t115);
t130 = t114 * t122;
t117 = sin(pkin(6));
t120 = sin(qJ(3));
t129 = t117 * t120;
t128 = t117 * t122;
t123 = cos(qJ(2));
t127 = t117 * t123;
t119 = cos(pkin(6));
t121 = sin(qJ(2));
t126 = t119 * t121;
t125 = t119 * t123;
t124 = t122 * t123;
t118 = cos(pkin(10));
t116 = sin(pkin(10));
t111 = t119 * t120 + t121 * t128;
t110 = t119 * t122 - t121 * t129;
t109 = -t116 * t126 + t118 * t123;
t108 = t116 * t125 + t118 * t121;
t107 = t116 * t123 + t118 * t126;
t106 = t116 * t121 - t118 * t125;
t105 = t109 * t122 + t116 * t129;
t104 = -t109 * t120 + t116 * t128;
t103 = t107 * t122 - t118 * t129;
t102 = -t107 * t120 - t118 * t128;
t1 = [0, -t108 * t130 + t109 * t113, t104 * t114, 0, -t105 * t113 + t108 * t114, 0; 0, -t106 * t130 + t107 * t113, t102 * t114, 0, -t103 * t113 + t106 * t114, 0; 0 (t113 * t121 + t114 * t124) * t117, t110 * t114, 0, -t111 * t113 - t114 * t127, 0; 0, t108 * t131 + t109 * t114, -t104 * t113, 0, -t105 * t114 - t108 * t113, 0; 0, t106 * t131 + t107 * t114, -t102 * t113, 0, -t103 * t114 - t106 * t113, 0; 0 (-t113 * t124 + t114 * t121) * t117, -t110 * t113, 0, -t111 * t114 + t113 * t127, 0; 0, -t108 * t120, t105, 0, 0, 0; 0, -t106 * t120, t103, 0, 0, 0; 0, t120 * t127, t111, 0, 0, 0;];
JR_rot  = t1;
