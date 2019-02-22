% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:11
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRP4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:11:25
% EndTime: 2019-02-22 11:11:25
% DurationCPUTime: 0.06s
% Computational Cost: add. (36->19), mult. (108->32), div. (0->0), fcn. (161->8), ass. (0->23)
t124 = sin(qJ(2));
t121 = sin(pkin(9));
t122 = cos(pkin(9));
t123 = sin(qJ(5));
t126 = cos(qJ(5));
t131 = t121 * t123 + t122 * t126;
t130 = t131 * t124;
t125 = sin(qJ(1));
t127 = cos(qJ(2));
t136 = t125 * t127;
t128 = cos(qJ(1));
t135 = t128 * t127;
t116 = t121 * t136 + t128 * t122;
t117 = -t128 * t121 + t122 * t136;
t134 = t116 * t126 - t117 * t123;
t133 = t116 * t123 + t117 * t126;
t132 = t121 * t126 - t122 * t123;
t129 = t132 * t124;
t119 = t125 * t121 + t122 * t135;
t118 = t121 * t135 - t125 * t122;
t115 = t118 * t123 + t119 * t126;
t114 = -t118 * t126 + t119 * t123;
t1 = [-t133, -t128 * t130, 0, 0, -t114, 0; t115, -t125 * t130, 0, 0, t134, 0; 0, t131 * t127, 0, 0, t129, 0; t125 * t124, -t135, 0, 0, 0, 0; -t128 * t124, -t136, 0, 0, 0, 0; 0, -t124, 0, 0, 0, 0; t134, t128 * t129, 0, 0, t115, 0; t114, t125 * t129, 0, 0, t133, 0; 0, -t132 * t127, 0, 0, t130, 0;];
JR_rot  = t1;
