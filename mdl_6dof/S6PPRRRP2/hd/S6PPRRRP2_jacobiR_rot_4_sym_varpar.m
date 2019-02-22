% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:25
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRP2_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobiR_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:25:25
% EndTime: 2019-02-22 09:25:25
% DurationCPUTime: 0.10s
% Computational Cost: add. (69->26), mult. (203->58), div. (0->0), fcn. (284->12), ass. (0->34)
t116 = sin(pkin(11));
t122 = cos(pkin(6));
t134 = t116 * t122;
t117 = sin(pkin(7));
t118 = sin(pkin(6));
t133 = t117 * t118;
t132 = t117 * t122;
t121 = cos(pkin(7));
t131 = t118 * t121;
t119 = cos(pkin(12));
t130 = t119 * t121;
t120 = cos(pkin(11));
t129 = t120 * t122;
t115 = sin(pkin(12));
t111 = -t116 * t115 + t119 * t129;
t128 = t111 * t121 - t120 * t133;
t113 = -t120 * t115 - t119 * t134;
t127 = t113 * t121 + t116 * t133;
t126 = cos(qJ(3));
t125 = cos(qJ(4));
t124 = sin(qJ(3));
t123 = sin(qJ(4));
t114 = -t115 * t134 + t120 * t119;
t112 = t115 * t129 + t116 * t119;
t110 = -t119 * t133 + t122 * t121;
t109 = -t113 * t117 + t116 * t131;
t108 = -t111 * t117 - t120 * t131;
t107 = t124 * t132 + (t115 * t126 + t124 * t130) * t118;
t106 = t126 * t132 + (-t115 * t124 + t126 * t130) * t118;
t105 = t114 * t126 + t127 * t124;
t104 = -t114 * t124 + t127 * t126;
t103 = t112 * t126 + t128 * t124;
t102 = -t112 * t124 + t128 * t126;
t1 = [0, 0, t104 * t125, -t105 * t123 + t109 * t125, 0, 0; 0, 0, t102 * t125, -t103 * t123 + t108 * t125, 0, 0; 0, 0, t106 * t125, -t107 * t123 + t110 * t125, 0, 0; 0, 0, -t104 * t123, -t105 * t125 - t109 * t123, 0, 0; 0, 0, -t102 * t123, -t103 * t125 - t108 * t123, 0, 0; 0, 0, -t106 * t123, -t107 * t125 - t110 * t123, 0, 0; 0, 0, t105, 0, 0, 0; 0, 0, t103, 0, 0, 0; 0, 0, t107, 0, 0, 0;];
JR_rot  = t1;
