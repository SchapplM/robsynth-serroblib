% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:21
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:21:06
% EndTime: 2019-02-22 12:21:07
% DurationCPUTime: 0.09s
% Computational Cost: add. (115->18), mult. (117->30), div. (0->0), fcn. (180->8), ass. (0->28)
t113 = sin(pkin(6));
t115 = sin(qJ(2));
t127 = t113 * t115;
t116 = sin(qJ(1));
t126 = t113 * t116;
t117 = cos(qJ(2));
t125 = t113 * t117;
t118 = cos(qJ(1));
t124 = t113 * t118;
t123 = t115 * t118;
t122 = t116 * t115;
t121 = t116 * t117;
t120 = t117 * t118;
t114 = cos(pkin(6));
t106 = t114 * t123 + t121;
t112 = qJ(3) + qJ(4) + pkin(12);
t110 = sin(t112);
t111 = cos(t112);
t100 = -t106 * t111 + t110 * t124;
t119 = t106 * t110 + t111 * t124;
t108 = -t114 * t122 + t120;
t107 = t114 * t121 + t123;
t105 = t114 * t120 - t122;
t104 = -t110 * t114 - t111 * t127;
t103 = -t110 * t127 + t111 * t114;
t102 = t108 * t111 + t110 * t126;
t101 = -t108 * t110 + t111 * t126;
t1 = [t100, -t107 * t111, t101, t101, 0, 0; t102, t105 * t111, -t119, -t119, 0, 0; 0, t111 * t125, t103, t103, 0, 0; t119, t107 * t110, -t102, -t102, 0, 0; t101, -t105 * t110, t100, t100, 0, 0; 0, -t110 * t125, t104, t104, 0, 0; t105, t108, 0, 0, 0, 0; t107, t106, 0, 0, 0, 0; 0, t127, 0, 0, 0, 0;];
JR_rot  = t1;
