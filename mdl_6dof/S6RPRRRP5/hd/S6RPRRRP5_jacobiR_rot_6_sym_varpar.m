% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:53
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRP5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:53:20
% EndTime: 2019-02-22 10:53:20
% DurationCPUTime: 0.05s
% Computational Cost: add. (78->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t116 = sin(qJ(5));
t117 = sin(qJ(1));
t127 = t117 * t116;
t118 = cos(qJ(5));
t126 = t117 * t118;
t119 = cos(qJ(1));
t125 = t119 * t116;
t124 = t119 * t118;
t115 = pkin(10) + qJ(3) + qJ(4);
t113 = sin(t115);
t123 = t113 * t127;
t122 = t113 * t126;
t121 = t113 * t125;
t120 = t113 * t124;
t114 = cos(t115);
t112 = t119 * t114;
t111 = t114 * t118;
t110 = t114 * t116;
t109 = t117 * t114;
t108 = t114 * t124 + t127;
t107 = t114 * t125 - t126;
t106 = t114 * t126 - t125;
t105 = -t114 * t127 - t124;
t1 = [-t106, 0, -t120, -t120, -t107, 0; t108, 0, -t122, -t122, t105, 0; 0, 0, t111, t111, -t113 * t116, 0; -t117 * t113, 0, t112, t112, 0, 0; t119 * t113, 0, t109, t109, 0, 0; 0, 0, t113, t113, 0, 0; t105, 0, -t121, -t121, t108, 0; t107, 0, -t123, -t123, t106, 0; 0, 0, t110, t110, t113 * t118, 0;];
JR_rot  = t1;
