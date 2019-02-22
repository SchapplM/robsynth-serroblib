% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:31
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:31:24
% EndTime: 2019-02-22 09:31:24
% DurationCPUTime: 0.07s
% Computational Cost: add. (90->28), mult. (163->65), div. (0->0), fcn. (238->10), ass. (0->31)
t117 = pkin(11) + qJ(4);
t115 = sin(t117);
t122 = sin(qJ(6));
t134 = t115 * t122;
t124 = cos(qJ(6));
t133 = t115 * t124;
t118 = sin(pkin(10));
t119 = sin(pkin(6));
t132 = t118 * t119;
t120 = cos(pkin(10));
t131 = t119 * t120;
t123 = sin(qJ(2));
t130 = t119 * t123;
t121 = cos(pkin(6));
t129 = t121 * t123;
t125 = cos(qJ(2));
t128 = t121 * t125;
t127 = t122 * t125;
t126 = t124 * t125;
t116 = cos(t117);
t113 = -t118 * t129 + t120 * t125;
t112 = t118 * t128 + t120 * t123;
t111 = t118 * t125 + t120 * t129;
t110 = t118 * t123 - t120 * t128;
t109 = t121 * t115 + t116 * t130;
t108 = t115 * t130 - t121 * t116;
t107 = t113 * t116 + t115 * t132;
t106 = t113 * t115 - t116 * t132;
t105 = t111 * t116 - t115 * t131;
t104 = t111 * t115 + t116 * t131;
t1 = [0, -t112 * t134 + t113 * t124, 0, t107 * t122, 0, t106 * t124 - t112 * t122; 0, -t110 * t134 + t111 * t124, 0, t105 * t122, 0, t104 * t124 - t110 * t122; 0 (t115 * t127 + t123 * t124) * t119, 0, t109 * t122, 0, t108 * t124 + t119 * t127; 0, -t112 * t133 - t113 * t122, 0, t107 * t124, 0, -t106 * t122 - t112 * t124; 0, -t110 * t133 - t111 * t122, 0, t105 * t124, 0, -t104 * t122 - t110 * t124; 0 (t115 * t126 - t122 * t123) * t119, 0, t109 * t124, 0, -t108 * t122 + t119 * t126; 0, -t112 * t116, 0, -t106, 0, 0; 0, -t110 * t116, 0, -t104, 0, 0; 0, t119 * t125 * t116, 0, -t108, 0, 0;];
JR_rot  = t1;
