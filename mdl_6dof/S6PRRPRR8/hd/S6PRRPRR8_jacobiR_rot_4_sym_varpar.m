% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR8_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiR_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:17
% EndTime: 2019-02-26 20:08:17
% DurationCPUTime: 0.06s
% Computational Cost: add. (38->20), mult. (118->51), div. (0->0), fcn. (169->10), ass. (0->27)
t116 = sin(pkin(7));
t117 = sin(pkin(6));
t136 = t116 * t117;
t120 = cos(pkin(6));
t135 = t116 * t120;
t119 = cos(pkin(7));
t121 = sin(qJ(3));
t134 = t119 * t121;
t123 = cos(qJ(3));
t133 = t119 * t123;
t122 = sin(qJ(2));
t132 = t120 * t122;
t124 = cos(qJ(2));
t131 = t120 * t124;
t130 = t121 * t122;
t129 = t121 * t124;
t128 = t122 * t123;
t127 = t123 * t124;
t115 = sin(pkin(12));
t118 = cos(pkin(12));
t110 = -t115 * t122 + t118 * t131;
t126 = -t110 * t119 + t118 * t136;
t112 = -t115 * t131 - t118 * t122;
t125 = t112 * t119 + t115 * t136;
t113 = -t115 * t132 + t118 * t124;
t111 = t115 * t124 + t118 * t132;
t1 = [0, t113 * t116, 0, 0, 0, 0; 0, t111 * t116, 0, 0, 0, 0; 0, t122 * t136, 0, 0, 0, 0; 0, -t112 * t123 + t113 * t134, t113 * t121 - t125 * t123, 0, 0, 0; 0, -t110 * t123 + t111 * t134, t111 * t121 + t126 * t123, 0, 0, 0; 0 (t119 * t130 - t127) * t117, -t123 * t135 + (-t119 * t127 + t130) * t117, 0, 0, 0; 0, t112 * t121 + t113 * t133, t113 * t123 + t125 * t121, 0, 0, 0; 0, t110 * t121 + t111 * t133, t111 * t123 - t126 * t121, 0, 0, 0; 0 (t119 * t128 + t129) * t117, t121 * t135 + (t119 * t129 + t128) * t117, 0, 0, 0;];
JR_rot  = t1;
