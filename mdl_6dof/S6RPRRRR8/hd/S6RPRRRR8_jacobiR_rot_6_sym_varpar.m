% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRRR8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:18:40
% EndTime: 2019-02-26 21:18:40
% DurationCPUTime: 0.08s
% Computational Cost: add. (102->23), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
t116 = qJ(3) + qJ(4);
t112 = sin(t116);
t115 = qJ(5) + qJ(6);
t113 = cos(t115);
t128 = t112 * t113;
t111 = sin(t115);
t114 = cos(t116);
t127 = t114 * t111;
t126 = t114 * t113;
t117 = sin(qJ(1));
t110 = t117 * t112;
t125 = t117 * t113;
t124 = t117 * t114;
t118 = cos(qJ(1));
t123 = t118 * t112;
t122 = t118 * t113;
t121 = t118 * t114;
t120 = t111 * t124;
t119 = t113 * t121;
t109 = t112 * t111;
t108 = t111 * t121;
t107 = t113 * t124;
t106 = -t111 * t117 + t112 * t122;
t105 = t111 * t123 + t125;
t104 = t111 * t118 + t112 * t125;
t103 = -t110 * t111 + t122;
t1 = [t106, 0, t107, t107, t103, t103; t104, 0, -t119, -t119, t105, t105; 0, 0, -t128, -t128, -t127, -t127; -t105, 0, -t120, -t120, -t104, -t104; t103, 0, t108, t108, t106, t106; 0, 0, t109, t109, -t126, -t126; -t121, 0, t110, t110, 0, 0; -t124, 0, -t123, -t123, 0, 0; 0, 0, t114, t114, 0, 0;];
JR_rot  = t1;
