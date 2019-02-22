% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:55
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR9_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobiR_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:55:00
% EndTime: 2019-02-22 11:55:00
% DurationCPUTime: 0.07s
% Computational Cost: add. (54->25), mult. (163->57), div. (0->0), fcn. (238->10), ass. (0->30)
t111 = sin(pkin(11));
t118 = cos(qJ(3));
t130 = t111 * t118;
t112 = sin(pkin(6));
t115 = sin(qJ(3));
t129 = t112 * t115;
t128 = t112 * t118;
t120 = cos(qJ(1));
t127 = t112 * t120;
t113 = cos(pkin(11));
t126 = t113 * t118;
t116 = sin(qJ(2));
t117 = sin(qJ(1));
t125 = t117 * t116;
t119 = cos(qJ(2));
t124 = t117 * t119;
t123 = t118 * t119;
t122 = t120 * t116;
t121 = t120 * t119;
t114 = cos(pkin(6));
t107 = t114 * t122 + t124;
t101 = -t107 * t115 - t118 * t127;
t102 = -t107 * t118 + t115 * t127;
t109 = -t114 * t125 + t121;
t108 = t114 * t124 + t122;
t106 = t114 * t121 - t125;
t105 = t114 * t118 - t116 * t129;
t104 = t109 * t118 + t117 * t129;
t103 = t109 * t115 - t117 * t128;
t1 = [t102 * t113 + t106 * t111, -t108 * t126 + t109 * t111, -t103 * t113, 0, 0, 0; t104 * t113 + t108 * t111, t106 * t126 + t107 * t111, t101 * t113, 0, 0, 0; 0 (t111 * t116 + t113 * t123) * t112, t105 * t113, 0, 0, 0; -t102 * t111 + t106 * t113, t108 * t130 + t109 * t113, t103 * t111, 0, 0, 0; -t104 * t111 + t108 * t113, -t106 * t130 + t107 * t113, -t101 * t111, 0, 0, 0; 0 (-t111 * t123 + t113 * t116) * t112, -t105 * t111, 0, 0, 0; t101, -t108 * t115, t104, 0, 0, 0; t103, t106 * t115, -t102, 0, 0, 0; 0, t119 * t129, t114 * t115 + t116 * t128, 0, 0, 0;];
JR_rot  = t1;
