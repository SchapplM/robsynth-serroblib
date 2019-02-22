% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:52
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:52:13
% EndTime: 2019-02-22 11:52:13
% DurationCPUTime: 0.07s
% Computational Cost: add. (93->26), mult. (163->58), div. (0->0), fcn. (238->10), ass. (0->31)
t121 = qJ(3) + pkin(11);
t120 = cos(t121);
t122 = sin(pkin(12));
t139 = t120 * t122;
t124 = cos(pkin(12));
t138 = t120 * t124;
t128 = cos(qJ(2));
t137 = t120 * t128;
t123 = sin(pkin(6));
t126 = sin(qJ(2));
t136 = t123 * t126;
t127 = sin(qJ(1));
t135 = t123 * t127;
t129 = cos(qJ(1));
t134 = t123 * t129;
t133 = t127 * t126;
t132 = t127 * t128;
t131 = t129 * t126;
t130 = t129 * t128;
t125 = cos(pkin(6));
t115 = t125 * t131 + t132;
t119 = sin(t121);
t109 = -t115 * t119 - t120 * t134;
t110 = -t115 * t120 + t119 * t134;
t117 = -t125 * t133 + t130;
t116 = t125 * t132 + t131;
t114 = t125 * t130 - t133;
t113 = -t119 * t136 + t125 * t120;
t112 = t117 * t120 + t119 * t135;
t111 = t117 * t119 - t120 * t135;
t1 = [t110 * t124 + t114 * t122, -t116 * t138 + t117 * t122, -t111 * t124, 0, 0, 0; t112 * t124 + t116 * t122, t114 * t138 + t115 * t122, t109 * t124, 0, 0, 0; 0 (t122 * t126 + t124 * t137) * t123, t113 * t124, 0, 0, 0; -t110 * t122 + t114 * t124, t116 * t139 + t117 * t124, t111 * t122, 0, 0, 0; -t112 * t122 + t116 * t124, -t114 * t139 + t115 * t124, -t109 * t122, 0, 0, 0; 0 (-t122 * t137 + t124 * t126) * t123, -t113 * t122, 0, 0, 0; t109, -t116 * t119, t112, 0, 0, 0; t111, t114 * t119, -t110, 0, 0, 0; 0, t123 * t128 * t119, t125 * t119 + t120 * t136, 0, 0, 0;];
JR_rot  = t1;
