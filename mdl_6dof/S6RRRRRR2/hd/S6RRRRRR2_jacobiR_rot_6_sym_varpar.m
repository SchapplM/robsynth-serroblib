% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:34
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:33:59
% EndTime: 2019-02-22 12:33:59
% DurationCPUTime: 0.08s
% Computational Cost: add. (162->24), mult. (78->20), div. (0->0), fcn. (135->6), ass. (0->27)
t121 = qJ(2) + qJ(3) + qJ(4);
t117 = sin(t121);
t122 = qJ(5) + qJ(6);
t119 = sin(t122);
t133 = t117 * t119;
t120 = cos(t122);
t132 = t117 * t120;
t118 = cos(t121);
t131 = t118 * t119;
t123 = sin(qJ(1));
t130 = t123 * t119;
t129 = t123 * t120;
t124 = cos(qJ(1));
t128 = t124 * t119;
t127 = t124 * t120;
t126 = t117 * t129;
t125 = t117 * t127;
t116 = t124 * t118;
t115 = t123 * t118;
t114 = t118 * t120;
t113 = t117 * t128;
t112 = t117 * t130;
t111 = t118 * t127 + t130;
t110 = -t118 * t128 + t129;
t109 = -t118 * t129 + t128;
t108 = t118 * t130 + t127;
t1 = [t109, -t125, -t125, -t125, t110, t110; t111, -t126, -t126, -t126, -t108, -t108; 0, t114, t114, t114, -t133, -t133; t108, t113, t113, t113, -t111, -t111; t110, t112, t112, t112, t109, t109; 0, -t131, -t131, -t131, -t132, -t132; -t123 * t117, t116, t116, t116, 0, 0; t124 * t117, t115, t115, t115, 0, 0; 0, t117, t117, t117, 0, 0;];
JR_rot  = t1;
