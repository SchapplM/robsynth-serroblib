% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR11_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:50
% EndTime: 2019-02-26 22:21:50
% DurationCPUTime: 0.09s
% Computational Cost: add. (91->32), mult. (265->47), div. (0->0), fcn. (382->10), ass. (0->37)
t120 = cos(pkin(6));
t123 = sin(qJ(2));
t128 = cos(qJ(1));
t136 = t128 * t123;
t124 = sin(qJ(1));
t127 = cos(qJ(2));
t137 = t124 * t127;
t113 = t120 * t136 + t137;
t122 = sin(qJ(3));
t126 = cos(qJ(3));
t119 = sin(pkin(6));
t139 = t119 * t128;
t104 = t113 * t122 + t126 * t139;
t107 = -t113 * t126 + t122 * t139;
t121 = sin(qJ(5));
t125 = cos(qJ(5));
t134 = t104 * t125 + t107 * t121;
t133 = t104 * t121 - t107 * t125;
t142 = t119 * t123;
t141 = t119 * t126;
t140 = t119 * t127;
t138 = t124 * t123;
t135 = t128 * t127;
t115 = -t120 * t138 + t135;
t108 = t115 * t122 - t124 * t141;
t109 = t124 * t119 * t122 + t115 * t126;
t102 = t108 * t125 - t109 * t121;
t103 = t108 * t121 + t109 * t125;
t110 = -t120 * t126 + t122 * t142;
t111 = t120 * t122 + t123 * t141;
t132 = t110 * t125 - t111 * t121;
t131 = t110 * t121 + t111 * t125;
t130 = t121 * t126 - t122 * t125;
t129 = t121 * t122 + t125 * t126;
t114 = -t120 * t137 - t136;
t112 = -t120 * t135 + t138;
t1 = [-t133, t129 * t114, -t102, 0, t102, 0; t103, -t129 * t112, -t134, 0, t134, 0; 0, t129 * t140, -t132, 0, t132, 0; -t134, -t130 * t114, t103, 0, -t103, 0; t102, t130 * t112, t133, 0, -t133, 0; 0, -t130 * t140, t131, 0, -t131, 0; t112, -t115, 0, 0, 0, 0; t114, -t113, 0, 0, 0, 0; 0, -t142, 0, 0, 0, 0;];
JR_rot  = t1;
