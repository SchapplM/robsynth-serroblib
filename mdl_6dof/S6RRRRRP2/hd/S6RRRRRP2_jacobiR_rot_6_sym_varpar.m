% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:40:10
% EndTime: 2019-02-26 22:40:10
% DurationCPUTime: 0.04s
% Computational Cost: add. (100->21), mult. (64->20), div. (0->0), fcn. (111->6), ass. (0->24)
t121 = sin(qJ(5));
t122 = sin(qJ(1));
t132 = t122 * t121;
t123 = cos(qJ(5));
t131 = t122 * t123;
t124 = cos(qJ(1));
t130 = t124 * t121;
t129 = t124 * t123;
t120 = qJ(2) + qJ(3) + qJ(4);
t118 = sin(t120);
t128 = t118 * t132;
t127 = t118 * t131;
t126 = t118 * t130;
t125 = t118 * t129;
t119 = cos(t120);
t117 = t124 * t119;
t116 = t119 * t123;
t115 = t119 * t121;
t114 = t122 * t119;
t113 = t119 * t129 + t132;
t112 = t119 * t130 - t131;
t111 = t119 * t131 - t130;
t110 = -t119 * t132 - t129;
t1 = [-t111, -t125, -t125, -t125, -t112, 0; t113, -t127, -t127, -t127, t110, 0; 0, t116, t116, t116, -t118 * t121, 0; -t122 * t118, t117, t117, t117, 0, 0; t124 * t118, t114, t114, t114, 0, 0; 0, t118, t118, t118, 0, 0; t110, -t126, -t126, -t126, t113, 0; t112, -t128, -t128, -t128, t111, 0; 0, t115, t115, t115, t118 * t123, 0;];
JR_rot  = t1;
