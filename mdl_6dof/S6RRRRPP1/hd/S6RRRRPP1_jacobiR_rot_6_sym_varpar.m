% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:53
% EndTime: 2019-02-26 22:24:53
% DurationCPUTime: 0.08s
% Computational Cost: add. (78->18), mult. (52->20), div. (0->0), fcn. (90->6), ass. (0->23)
t118 = qJ(2) + qJ(3);
t115 = sin(t118);
t119 = sin(qJ(1));
t126 = t119 * t115;
t116 = cos(t118);
t111 = t119 * t116;
t120 = cos(qJ(1));
t125 = t120 * t115;
t112 = t120 * t116;
t117 = qJ(4) + pkin(10);
t113 = sin(t117);
t124 = t113 * t126;
t114 = cos(t117);
t123 = t114 * t126;
t122 = t113 * t125;
t121 = t114 * t125;
t110 = t116 * t114;
t109 = t116 * t113;
t108 = t114 * t112 + t119 * t113;
t107 = t113 * t112 - t119 * t114;
t106 = t114 * t111 - t120 * t113;
t105 = -t113 * t111 - t120 * t114;
t1 = [-t106, -t121, -t121, -t107, 0, 0; t108, -t123, -t123, t105, 0, 0; 0, t110, t110, -t115 * t113, 0, 0; -t126, t112, t112, 0, 0, 0; t125, t111, t111, 0, 0, 0; 0, t115, t115, 0, 0, 0; t105, -t122, -t122, t108, 0, 0; t107, -t124, -t124, t106, 0, 0; 0, t109, t109, t115 * t114, 0, 0;];
JR_rot  = t1;
