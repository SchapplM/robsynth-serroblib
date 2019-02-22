% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:25
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR5_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobiR_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:25:14
% EndTime: 2019-02-22 11:25:14
% DurationCPUTime: 0.06s
% Computational Cost: add. (66->18), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->27)
t115 = sin(pkin(6));
t120 = sin(qJ(1));
t128 = t115 * t120;
t123 = cos(qJ(1));
t127 = t115 * t123;
t117 = cos(pkin(6));
t114 = sin(pkin(11));
t116 = cos(pkin(11));
t119 = sin(qJ(2));
t122 = cos(qJ(2));
t125 = t122 * t114 + t119 * t116;
t109 = t125 * t117;
t110 = t119 * t114 - t122 * t116;
t103 = t123 * t109 - t120 * t110;
t118 = sin(qJ(4));
t121 = cos(qJ(4));
t126 = -t103 * t121 + t118 * t127;
t105 = -t120 * t109 - t123 * t110;
t124 = t103 * t118 + t121 * t127;
t108 = t110 * t117;
t107 = t125 * t115;
t106 = t110 * t115;
t104 = t120 * t108 - t123 * t125;
t102 = -t123 * t108 - t120 * t125;
t101 = t105 * t121 + t118 * t128;
t100 = -t105 * t118 + t121 * t128;
t1 = [t126, t104 * t121, 0, t100, 0, 0; t101, t102 * t121, 0, -t124, 0, 0; 0, -t106 * t121, 0, -t107 * t118 + t117 * t121, 0, 0; t124, -t104 * t118, 0, -t101, 0, 0; t100, -t102 * t118, 0, t126, 0, 0; 0, t106 * t118, 0, -t107 * t121 - t117 * t118, 0, 0; t102, t105, 0, 0, 0, 0; -t104, t103, 0, 0, 0, 0; 0, t107, 0, 0, 0, 0;];
JR_rot  = t1;
