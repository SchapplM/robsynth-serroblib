% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:49:26
% EndTime: 2019-02-26 22:49:26
% DurationCPUTime: 0.11s
% Computational Cost: add. (149->20), mult. (147->30), div. (0->0), fcn. (226->8), ass. (0->28)
t113 = sin(pkin(6));
t115 = sin(qJ(2));
t127 = t113 * t115;
t116 = sin(qJ(1));
t126 = t113 * t116;
t117 = cos(qJ(2));
t125 = t113 * t117;
t118 = cos(qJ(1));
t124 = t113 * t118;
t123 = t116 * t115;
t122 = t116 * t117;
t121 = t118 * t115;
t120 = t118 * t117;
t114 = cos(pkin(6));
t106 = t114 * t121 + t122;
t112 = qJ(3) + qJ(4) + qJ(5);
t110 = sin(t112);
t111 = cos(t112);
t100 = -t106 * t111 + t110 * t124;
t119 = t106 * t110 + t111 * t124;
t108 = -t114 * t123 + t120;
t107 = t114 * t122 + t121;
t105 = t114 * t120 - t123;
t104 = -t114 * t110 - t111 * t127;
t103 = -t110 * t127 + t114 * t111;
t102 = t108 * t111 + t110 * t126;
t101 = -t108 * t110 + t111 * t126;
t1 = [t100, -t107 * t111, t101, t101, t101, 0; t102, t105 * t111, -t119, -t119, -t119, 0; 0, t111 * t125, t103, t103, t103, 0; t119, t107 * t110, -t102, -t102, -t102, 0; t101, -t105 * t110, t100, t100, t100, 0; 0, -t110 * t125, t104, t104, t104, 0; t105, t108, 0, 0, 0, 0; t107, t106, 0, 0, 0, 0; 0, t127, 0, 0, 0, 0;];
JR_rot  = t1;
