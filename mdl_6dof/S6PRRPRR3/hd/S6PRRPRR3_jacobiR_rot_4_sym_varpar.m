% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiR_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:24
% EndTime: 2019-02-26 20:05:24
% DurationCPUTime: 0.10s
% Computational Cost: add. (68->24), mult. (196->58), div. (0->0), fcn. (277->12), ass. (0->27)
t109 = sin(pkin(12));
t111 = sin(pkin(6));
t124 = t109 * t111;
t113 = cos(pkin(12));
t123 = t111 * t113;
t115 = cos(pkin(6));
t117 = sin(qJ(2));
t122 = t115 * t117;
t119 = cos(qJ(2));
t121 = t115 * t119;
t108 = sin(pkin(13));
t112 = cos(pkin(13));
t116 = sin(qJ(3));
t118 = cos(qJ(3));
t120 = t118 * t108 + t116 * t112;
t104 = t116 * t108 - t118 * t112;
t114 = cos(pkin(7));
t110 = sin(pkin(7));
t103 = -t109 * t122 + t113 * t119;
t102 = -t109 * t121 - t113 * t117;
t101 = t109 * t119 + t113 * t122;
t100 = -t109 * t117 + t113 * t121;
t99 = t120 * t114;
t98 = t104 * t114;
t97 = t120 * t110;
t96 = t104 * t110;
t1 = [0, -t102 * t104 - t103 * t99, -t102 * t98 - t103 * t120 - t96 * t124, 0, 0, 0; 0, -t100 * t104 - t101 * t99, -t100 * t98 - t101 * t120 + t96 * t123, 0, 0, 0; 0 (-t104 * t119 - t117 * t99) * t111, -t115 * t96 + (-t117 * t120 - t119 * t98) * t111, 0, 0, 0; 0, -t102 * t120 + t103 * t98, -t102 * t99 + t103 * t104 - t97 * t124, 0, 0, 0; 0, -t100 * t120 + t101 * t98, -t100 * t99 + t101 * t104 + t97 * t123, 0, 0, 0; 0 (t117 * t98 - t119 * t120) * t111, -t115 * t97 + (t104 * t117 - t119 * t99) * t111, 0, 0, 0; 0, t103 * t110, 0, 0, 0, 0; 0, t101 * t110, 0, 0, 0, 0; 0, t111 * t117 * t110, 0, 0, 0, 0;];
JR_rot  = t1;
