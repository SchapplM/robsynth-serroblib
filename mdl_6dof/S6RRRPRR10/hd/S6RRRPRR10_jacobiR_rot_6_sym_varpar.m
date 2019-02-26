% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:23
% EndTime: 2019-02-26 22:21:23
% DurationCPUTime: 0.11s
% Computational Cost: add. (120->30), mult. (188->32), div. (0->0), fcn. (281->8), ass. (0->24)
t116 = sin(qJ(2));
t114 = qJ(5) + qJ(6);
t112 = sin(t114);
t113 = cos(t114);
t115 = sin(qJ(3));
t118 = cos(qJ(3));
t124 = t112 * t118 - t113 * t115;
t133 = t116 * t124;
t123 = t112 * t115 + t113 * t118;
t132 = t123 * t116;
t120 = cos(qJ(1));
t117 = sin(qJ(1));
t119 = cos(qJ(2));
t128 = t117 * t119;
t105 = t115 * t128 + t120 * t118;
t106 = -t120 * t115 + t118 * t128;
t98 = t105 * t113 - t106 * t112;
t127 = t120 * t119;
t125 = t105 * t112 + t106 * t113;
t107 = t115 * t127 - t117 * t118;
t108 = t117 * t115 + t118 * t127;
t100 = t107 * t113 - t108 * t112;
t101 = t107 * t112 + t108 * t113;
t1 = [-t125, -t120 * t132, -t100, 0, t100, t100; t101, -t117 * t132, -t98, 0, t98, t98; 0, t123 * t119, t133, 0, -t133, -t133; -t98, t120 * t133, t101, 0, -t101, -t101; t100, t117 * t133, t125, 0, -t125, -t125; 0, -t124 * t119, t132, 0, -t132, -t132; t117 * t116, -t127, 0, 0, 0, 0; -t120 * t116, -t128, 0, 0, 0, 0; 0, -t116, 0, 0, 0, 0;];
JR_rot  = t1;
