% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:20
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR5_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:19:54
% EndTime: 2019-02-22 12:19:54
% DurationCPUTime: 0.04s
% Computational Cost: add. (48->17), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->24)
t115 = sin(qJ(4));
t116 = sin(qJ(1));
t126 = t116 * t115;
t117 = cos(qJ(4));
t125 = t116 * t117;
t118 = cos(qJ(1));
t124 = t118 * t115;
t123 = t118 * t117;
t114 = qJ(2) + qJ(3);
t112 = sin(t114);
t122 = t112 * t126;
t121 = t112 * t125;
t120 = t112 * t124;
t119 = t112 * t123;
t113 = cos(t114);
t111 = t118 * t113;
t110 = t113 * t117;
t109 = t113 * t115;
t108 = t116 * t113;
t107 = t113 * t123 + t126;
t106 = t113 * t124 - t125;
t105 = t113 * t125 - t124;
t104 = -t113 * t126 - t123;
t1 = [-t105, -t119, -t119, -t106, 0, 0; t107, -t121, -t121, t104, 0, 0; 0, t110, t110, -t112 * t115, 0, 0; -t116 * t112, t111, t111, 0, 0, 0; t118 * t112, t108, t108, 0, 0, 0; 0, t112, t112, 0, 0, 0; t104, -t120, -t120, t107, 0, 0; t106, -t122, -t122, t105, 0, 0; 0, t109, t109, t112 * t117, 0, 0;];
JR_rot  = t1;
