% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:48:18
% EndTime: 2019-02-26 22:48:18
% DurationCPUTime: 0.05s
% Computational Cost: add. (171->25), mult. (80->20), div. (0->0), fcn. (138->6), ass. (0->25)
t120 = qJ(4) + qJ(5) + qJ(6);
t116 = sin(t120);
t121 = qJ(2) + qJ(3);
t118 = sin(t121);
t130 = t118 * t116;
t117 = cos(t120);
t129 = t118 * t117;
t119 = cos(t121);
t128 = t119 * t116;
t122 = sin(qJ(1));
t127 = t122 * t118;
t114 = t122 * t119;
t123 = cos(qJ(1));
t126 = t123 * t118;
t115 = t123 * t119;
t125 = t117 * t127;
t124 = t117 * t126;
t113 = t119 * t117;
t112 = t116 * t126;
t111 = t116 * t127;
t110 = t117 * t115 + t122 * t116;
t109 = -t116 * t115 + t122 * t117;
t108 = -t117 * t114 + t123 * t116;
t107 = t116 * t114 + t123 * t117;
t1 = [t108, -t124, -t124, t109, t109, t109; t110, -t125, -t125, -t107, -t107, -t107; 0, t113, t113, -t130, -t130, -t130; t107, t112, t112, -t110, -t110, -t110; t109, t111, t111, t108, t108, t108; 0, -t128, -t128, -t129, -t129, -t129; -t127, t115, t115, 0, 0, 0; t126, t114, t114, 0, 0, 0; 0, t118, t118, 0, 0, 0;];
JR_rot  = t1;
