% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:46
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR12_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:46:32
% EndTime: 2019-02-22 11:46:32
% DurationCPUTime: 0.12s
% Computational Cost: add. (74->17), mult. (117->30), div. (0->0), fcn. (180->8), ass. (0->28)
t112 = sin(pkin(6));
t114 = sin(qJ(2));
t125 = t112 * t114;
t115 = sin(qJ(1));
t124 = t112 * t115;
t116 = cos(qJ(2));
t123 = t112 * t116;
t117 = cos(qJ(1));
t122 = t112 * t117;
t121 = t115 * t114;
t120 = t115 * t116;
t119 = t117 * t114;
t118 = t117 * t116;
t113 = cos(pkin(6));
t103 = -t113 * t118 + t121;
t111 = qJ(4) + qJ(5);
t109 = sin(t111);
t110 = cos(t111);
t100 = -t103 * t109 + t110 * t122;
t99 = t103 * t110 + t109 * t122;
t106 = -t113 * t121 + t118;
t105 = t113 * t120 + t119;
t104 = t113 * t119 + t120;
t102 = t109 * t123 - t113 * t110;
t101 = -t113 * t109 - t110 * t123;
t98 = t105 * t109 + t110 * t124;
t97 = t105 * t110 - t109 * t124;
t1 = [t100, t106 * t109, 0, t97, t97, 0; t98, t104 * t109, 0, t99, t99, 0; 0, t109 * t125, 0, t101, t101, 0; -t99, t106 * t110, 0, -t98, -t98, 0; t97, t104 * t110, 0, t100, t100, 0; 0, t110 * t125, 0, t102, t102, 0; -t104, -t105, 0, 0, 0, 0; t106, -t103, 0, 0, 0, 0; 0, t123, 0, 0, 0, 0;];
JR_rot  = t1;
