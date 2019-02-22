% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:37
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:37:39
% EndTime: 2019-02-22 09:37:39
% DurationCPUTime: 0.11s
% Computational Cost: add. (94->16), mult. (178->36), div. (0->0), fcn. (260->10), ass. (0->30)
t120 = sin(pkin(11));
t121 = sin(pkin(6));
t129 = t120 * t121;
t123 = cos(pkin(11));
t128 = t121 * t123;
t124 = cos(pkin(6));
t119 = sin(pkin(12));
t122 = cos(pkin(12));
t125 = sin(qJ(2));
t126 = cos(qJ(2));
t127 = t119 * t126 + t125 * t122;
t112 = t127 * t124;
t113 = t119 * t125 - t126 * t122;
t104 = t112 * t123 - t113 * t120;
t106 = -t112 * t120 - t113 * t123;
t118 = qJ(4) + qJ(5);
t117 = cos(t118);
t116 = sin(t118);
t111 = t113 * t124;
t110 = t127 * t121;
t109 = t113 * t121;
t108 = -t110 * t117 - t116 * t124;
t107 = -t110 * t116 + t117 * t124;
t105 = t111 * t120 - t123 * t127;
t103 = -t111 * t123 - t120 * t127;
t102 = -t106 * t117 - t116 * t129;
t101 = -t106 * t116 + t117 * t129;
t100 = -t104 * t117 + t116 * t128;
t99 = -t104 * t116 - t117 * t128;
t1 = [0, t105 * t117, 0, t101, t101, 0; 0, t103 * t117, 0, t99, t99, 0; 0, -t109 * t117, 0, t107, t107, 0; 0, -t105 * t116, 0, t102, t102, 0; 0, -t103 * t116, 0, t100, t100, 0; 0, t109 * t116, 0, t108, t108, 0; 0, t106, 0, 0, 0, 0; 0, t104, 0, 0, 0, 0; 0, t110, 0, 0, 0, 0;];
JR_rot  = t1;
