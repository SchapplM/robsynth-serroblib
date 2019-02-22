% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:32
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP12_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobiR_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobiR_rot_3_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:32:41
% EndTime: 2019-02-22 12:32:42
% DurationCPUTime: 0.11s
% Computational Cost: add. (58->25), mult. (178->58), div. (0->0), fcn. (255->10), ass. (0->32)
t115 = sin(qJ(2));
t116 = sin(qJ(1));
t118 = cos(qJ(2));
t119 = cos(qJ(1));
t136 = cos(pkin(6));
t124 = t119 * t136;
t104 = t116 * t115 - t118 * t124;
t105 = t115 * t124 + t116 * t118;
t113 = cos(pkin(7));
t114 = sin(qJ(3));
t117 = cos(qJ(3));
t111 = sin(pkin(7));
t112 = sin(pkin(6));
t133 = t112 * t119;
t126 = t111 * t133;
t137 = t117 * (t104 * t113 + t126) + t105 * t114;
t134 = t112 * t116;
t132 = t113 * t114;
t131 = t113 * t117;
t130 = t114 * t115;
t129 = t114 * t118;
t128 = t115 * t117;
t127 = t117 * t118;
t125 = t116 * t136;
t123 = t136 * t111;
t106 = -t119 * t115 - t118 * t125;
t121 = t106 * t113 + t111 * t134;
t120 = t104 * t132 - t105 * t117 + t114 * t126;
t107 = -t115 * t125 + t119 * t118;
t103 = t107 * t117 + t114 * t121;
t102 = -t107 * t114 + t117 * t121;
t1 = [t120, t106 * t117 - t107 * t132, t102, 0, 0, 0; t103, -t104 * t117 - t105 * t132, -t137, 0, 0, 0; 0 (-t113 * t130 + t127) * t112, t117 * t123 + (t113 * t127 - t130) * t112, 0, 0, 0; t137, -t106 * t114 - t107 * t131, -t103, 0, 0, 0; t102, t104 * t114 - t105 * t131, t120, 0, 0, 0; 0 (-t113 * t128 - t129) * t112, -t114 * t123 + (-t113 * t129 - t128) * t112, 0, 0, 0; -t104 * t111 + t113 * t133, t107 * t111, 0, 0, 0, 0; -t106 * t111 + t113 * t134, t105 * t111, 0, 0, 0, 0; 0, t112 * t115 * t111, 0, 0, 0, 0;];
JR_rot  = t1;
