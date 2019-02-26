% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:20
% EndTime: 2019-02-26 22:41:20
% DurationCPUTime: 0.04s
% Computational Cost: add. (99->21), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
t115 = qJ(4) + qJ(5);
t111 = sin(t115);
t116 = qJ(2) + qJ(3);
t112 = sin(t116);
t127 = t112 * t111;
t113 = cos(t115);
t126 = t112 * t113;
t114 = cos(t116);
t125 = t114 * t111;
t117 = sin(qJ(1));
t124 = t117 * t112;
t123 = t117 * t113;
t109 = t117 * t114;
t118 = cos(qJ(1));
t122 = t118 * t112;
t121 = t118 * t113;
t110 = t118 * t114;
t120 = t112 * t123;
t119 = t112 * t121;
t108 = t114 * t113;
t107 = t111 * t122;
t106 = t111 * t124;
t105 = t113 * t110 + t117 * t111;
t104 = -t111 * t110 + t123;
t103 = -t113 * t109 + t118 * t111;
t102 = t111 * t109 + t121;
t1 = [t103, -t119, -t119, t104, t104, 0; t105, -t120, -t120, -t102, -t102, 0; 0, t108, t108, -t127, -t127, 0; t102, t107, t107, -t105, -t105, 0; t104, t106, t106, t103, t103, 0; 0, -t125, -t125, -t126, -t126, 0; -t124, t110, t110, 0, 0, 0; t122, t109, t109, 0, 0, 0; 0, t112, t112, 0, 0, 0;];
JR_rot  = t1;
