% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPP1_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobiR_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:48
% EndTime: 2019-02-26 22:02:48
% DurationCPUTime: 0.15s
% Computational Cost: add. (53->27), mult. (183->71), div. (0->0), fcn. (258->10), ass. (0->35)
t100 = sin(pkin(10));
t103 = cos(pkin(6));
t129 = t100 * t103;
t101 = sin(pkin(6));
t105 = sin(qJ(2));
t128 = t101 * t105;
t102 = cos(pkin(10));
t127 = t102 * t103;
t107 = cos(qJ(3));
t126 = t103 * t107;
t104 = sin(qJ(3));
t108 = cos(qJ(2));
t125 = t104 * t108;
t124 = t105 * t103;
t106 = sin(qJ(1));
t123 = t105 * t106;
t122 = t105 * t107;
t109 = cos(qJ(1));
t121 = t105 * t109;
t120 = t106 * t104;
t119 = t107 * t108;
t118 = t109 * t104;
t117 = t109 * t107;
t96 = t108 * t120 + t117;
t116 = -t101 * t123 + t103 * t96;
t98 = t106 * t107 - t108 * t118;
t115 = t101 * t121 + t103 * t98;
t114 = -t103 * t125 + t128;
t113 = t101 * t108 + t104 * t124;
t112 = t103 * t108 - t104 * t128;
t111 = t100 * t122 + t113 * t102;
t110 = t113 * t100 - t102 * t122;
t99 = t108 * t117 + t120;
t97 = -t106 * t119 + t118;
t1 = [t116 * t100 + t97 * t102, t110 * t109, t98 * t102 - t99 * t129, 0, 0, 0; t115 * t100 + t99 * t102, t110 * t106, -t96 * t102 + t97 * t129, 0, 0, 0; 0, t114 * t100 + t102 * t119 (-t100 * t126 - t102 * t104) * t105, 0, 0, 0; -t97 * t100 + t116 * t102, t111 * t109, -t98 * t100 - t99 * t127, 0, 0, 0; -t99 * t100 + t115 * t102, t111 * t106, t96 * t100 + t97 * t127, 0, 0, 0; 0, -t100 * t119 + t114 * t102 (t100 * t104 - t102 * t126) * t105, 0, 0, 0; -t96 * t101 - t103 * t123, t112 * t109, t99 * t101, 0, 0, 0; -t98 * t101 + t103 * t121, t112 * t106, -t97 * t101, 0, 0, 0; 0, t101 * t125 + t124, t101 * t122, 0, 0, 0;];
JR_rot  = t1;
