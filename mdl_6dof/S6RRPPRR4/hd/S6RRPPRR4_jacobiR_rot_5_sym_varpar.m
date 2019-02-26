% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:11
% EndTime: 2019-02-26 21:30:11
% DurationCPUTime: 0.06s
% Computational Cost: add. (63->16), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->27)
t119 = sin(pkin(6));
t124 = sin(qJ(1));
t133 = t119 * t124;
t127 = cos(qJ(1));
t132 = t119 * t127;
t118 = sin(pkin(11));
t120 = cos(pkin(11));
t123 = sin(qJ(2));
t126 = cos(qJ(2));
t114 = t123 * t118 - t126 * t120;
t121 = cos(pkin(6));
t128 = t114 * t121;
t130 = t126 * t118 + t123 * t120;
t107 = -t124 * t130 - t127 * t128;
t122 = sin(qJ(5));
t125 = cos(qJ(5));
t131 = t107 * t122 + t125 * t132;
t113 = t130 * t121;
t106 = t127 * t113 - t124 * t114;
t108 = -t124 * t113 - t127 * t114;
t129 = -t107 * t125 + t122 * t132;
t112 = t130 * t119;
t111 = t114 * t119;
t109 = t124 * t128 - t127 * t130;
t105 = -t109 * t122 + t125 * t133;
t104 = -t109 * t125 - t122 * t133;
t1 = [t131, t108 * t122, 0, 0, t104, 0; t105, t106 * t122, 0, 0, t129, 0; 0, t112 * t122, 0, 0, t111 * t125 - t121 * t122, 0; -t129, t108 * t125, 0, 0, -t105, 0; t104, t106 * t125, 0, 0, t131, 0; 0, t112 * t125, 0, 0, -t111 * t122 - t121 * t125, 0; -t106, t109, 0, 0, 0, 0; t108, t107, 0, 0, 0, 0; 0, -t111, 0, 0, 0, 0;];
JR_rot  = t1;
