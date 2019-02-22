% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:05
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:05:20
% EndTime: 2019-02-22 12:05:20
% DurationCPUTime: 0.04s
% Computational Cost: add. (139->21), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->25)
t115 = pkin(11) + qJ(5) + qJ(6);
t113 = sin(t115);
t118 = qJ(2) + qJ(3);
t116 = sin(t118);
t127 = t116 * t113;
t114 = cos(t115);
t126 = t116 * t114;
t117 = cos(t118);
t125 = t117 * t113;
t119 = sin(qJ(1));
t124 = t119 * t116;
t111 = t119 * t117;
t120 = cos(qJ(1));
t123 = t120 * t116;
t112 = t120 * t117;
t122 = t114 * t124;
t121 = t114 * t123;
t110 = t117 * t114;
t109 = t113 * t123;
t108 = t113 * t124;
t107 = t114 * t112 + t119 * t113;
t106 = -t113 * t112 + t119 * t114;
t105 = -t114 * t111 + t120 * t113;
t104 = t113 * t111 + t120 * t114;
t1 = [t105, -t121, -t121, 0, t106, t106; t107, -t122, -t122, 0, -t104, -t104; 0, t110, t110, 0, -t127, -t127; t104, t109, t109, 0, -t107, -t107; t106, t108, t108, 0, t105, t105; 0, -t125, -t125, 0, -t126, -t126; -t124, t112, t112, 0, 0, 0; t123, t111, t111, 0, 0, 0; 0, t116, t116, 0, 0, 0;];
JR_rot  = t1;
