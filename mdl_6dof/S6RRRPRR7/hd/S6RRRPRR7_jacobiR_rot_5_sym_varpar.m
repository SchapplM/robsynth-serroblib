% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:07
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:07:04
% EndTime: 2019-02-22 12:07:04
% DurationCPUTime: 0.06s
% Computational Cost: add. (115->18), mult. (117->30), div. (0->0), fcn. (180->8), ass. (0->28)
t111 = sin(pkin(6));
t113 = sin(qJ(2));
t125 = t111 * t113;
t114 = sin(qJ(1));
t124 = t111 * t114;
t115 = cos(qJ(2));
t123 = t111 * t115;
t116 = cos(qJ(1));
t122 = t111 * t116;
t121 = t114 * t113;
t120 = t114 * t115;
t119 = t116 * t113;
t118 = t116 * t115;
t112 = cos(pkin(6));
t104 = t112 * t119 + t120;
t110 = qJ(3) + pkin(12) + qJ(5);
t108 = sin(t110);
t109 = cos(t110);
t98 = -t104 * t109 + t108 * t122;
t117 = t104 * t108 + t109 * t122;
t106 = -t112 * t121 + t118;
t105 = t112 * t120 + t119;
t103 = t112 * t118 - t121;
t102 = -t112 * t108 - t109 * t125;
t101 = -t108 * t125 + t112 * t109;
t100 = t106 * t109 + t108 * t124;
t99 = -t106 * t108 + t109 * t124;
t1 = [t98, -t105 * t109, t99, 0, t99, 0; t100, t103 * t109, -t117, 0, -t117, 0; 0, t109 * t123, t101, 0, t101, 0; t117, t105 * t108, -t100, 0, -t100, 0; t99, -t103 * t108, t98, 0, t98, 0; 0, -t108 * t123, t102, 0, t102, 0; t103, t106, 0, 0, 0, 0; t105, t104, 0, 0, 0, 0; 0, t125, 0, 0, 0, 0;];
JR_rot  = t1;
