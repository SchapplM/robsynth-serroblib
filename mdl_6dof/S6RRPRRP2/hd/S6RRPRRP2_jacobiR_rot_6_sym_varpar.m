% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:46:35
% EndTime: 2019-02-26 21:46:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (78->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t118 = sin(qJ(5));
t119 = sin(qJ(1));
t129 = t119 * t118;
t120 = cos(qJ(5));
t128 = t119 * t120;
t121 = cos(qJ(1));
t127 = t121 * t118;
t126 = t121 * t120;
t117 = qJ(2) + pkin(10) + qJ(4);
t115 = sin(t117);
t125 = t115 * t129;
t124 = t115 * t128;
t123 = t115 * t127;
t122 = t115 * t126;
t116 = cos(t117);
t114 = t121 * t116;
t113 = t116 * t120;
t112 = t116 * t118;
t111 = t119 * t116;
t110 = t116 * t126 + t129;
t109 = t116 * t127 - t128;
t108 = t116 * t128 - t127;
t107 = -t116 * t129 - t126;
t1 = [-t108, -t122, 0, -t122, -t109, 0; t110, -t124, 0, -t124, t107, 0; 0, t113, 0, t113, -t115 * t118, 0; -t119 * t115, t114, 0, t114, 0, 0; t121 * t115, t111, 0, t111, 0, 0; 0, t115, 0, t115, 0, 0; t107, -t123, 0, -t123, t110, 0; t109, -t125, 0, -t125, t108, 0; 0, t112, 0, t112, t115 * t120, 0;];
JR_rot  = t1;
