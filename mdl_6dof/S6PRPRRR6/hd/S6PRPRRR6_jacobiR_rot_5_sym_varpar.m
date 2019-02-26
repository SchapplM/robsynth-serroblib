% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:51
% EndTime: 2019-02-26 19:56:51
% DurationCPUTime: 0.07s
% Computational Cost: add. (57->28), mult. (163->64), div. (0->0), fcn. (238->10), ass. (0->30)
t108 = sin(pkin(6));
t112 = sin(qJ(4));
t125 = t108 * t112;
t115 = cos(qJ(4));
t124 = t108 * t115;
t116 = cos(qJ(2));
t123 = t108 * t116;
t110 = cos(pkin(6));
t113 = sin(qJ(2));
t122 = t110 * t113;
t121 = t110 * t116;
t111 = sin(qJ(5));
t120 = t111 * t112;
t119 = t111 * t113;
t114 = cos(qJ(5));
t118 = t112 * t114;
t117 = t113 * t114;
t109 = cos(pkin(11));
t107 = sin(pkin(11));
t105 = t110 * t115 - t112 * t123;
t104 = -t110 * t112 - t115 * t123;
t103 = -t107 * t122 + t109 * t116;
t102 = t107 * t121 + t109 * t113;
t101 = t107 * t116 + t109 * t122;
t100 = t107 * t113 - t109 * t121;
t99 = t100 * t112 - t109 * t124;
t98 = t100 * t115 + t109 * t125;
t97 = t102 * t112 + t107 * t124;
t96 = t102 * t115 - t107 * t125;
t1 = [0, -t102 * t111 + t103 * t118, 0, t96 * t114, t103 * t114 - t97 * t111, 0; 0, -t100 * t111 + t101 * t118, 0, t98 * t114, t101 * t114 - t99 * t111, 0; 0 (t111 * t116 + t112 * t117) * t108, 0, t104 * t114, -t105 * t111 + t108 * t117, 0; 0, -t102 * t114 - t103 * t120, 0, -t96 * t111, -t103 * t111 - t97 * t114, 0; 0, -t100 * t114 - t101 * t120, 0, -t98 * t111, -t101 * t111 - t99 * t114, 0; 0 (-t112 * t119 + t114 * t116) * t108, 0, -t104 * t111, -t105 * t114 - t108 * t119, 0; 0, -t103 * t115, 0, t97, 0, 0; 0, -t101 * t115, 0, t99, 0, 0; 0, -t113 * t124, 0, t105, 0, 0;];
JR_rot  = t1;
