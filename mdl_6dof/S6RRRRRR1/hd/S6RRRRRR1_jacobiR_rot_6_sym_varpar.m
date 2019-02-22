% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR1
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
% Datum: 2019-02-22 12:33
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:33:30
% EndTime: 2019-02-22 12:33:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (167->22), mult. (76->20), div. (0->0), fcn. (132->6), ass. (0->24)
t112 = qJ(2) + qJ(3) + qJ(4) + qJ(5);
t111 = cos(t112);
t113 = sin(qJ(6));
t123 = t111 * t113;
t114 = sin(qJ(1));
t122 = t114 * t113;
t115 = cos(qJ(6));
t121 = t114 * t115;
t116 = cos(qJ(1));
t120 = t116 * t113;
t119 = t116 * t115;
t110 = sin(t112);
t118 = t110 * t121;
t117 = t110 * t119;
t109 = t116 * t111;
t108 = t111 * t115;
t107 = t114 * t111;
t106 = t110 * t120;
t105 = t110 * t122;
t104 = t111 * t119 + t122;
t103 = -t111 * t120 + t121;
t102 = -t111 * t121 + t120;
t101 = t111 * t122 + t119;
t1 = [t102, -t117, -t117, -t117, -t117, t103; t104, -t118, -t118, -t118, -t118, -t101; 0, t108, t108, t108, t108, -t110 * t113; t101, t106, t106, t106, t106, -t104; t103, t105, t105, t105, t105, t102; 0, -t123, -t123, -t123, -t123, -t110 * t115; -t114 * t110, t109, t109, t109, t109, 0; t116 * t110, t107, t107, t107, t107, 0; 0, t110, t110, t110, t110, 0;];
JR_rot  = t1;
