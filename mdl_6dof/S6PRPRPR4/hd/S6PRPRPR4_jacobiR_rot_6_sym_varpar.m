% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:48:15
% EndTime: 2019-02-26 19:48:15
% DurationCPUTime: 0.07s
% Computational Cost: add. (123->29), mult. (163->64), div. (0->0), fcn. (238->10), ass. (0->32)
t123 = pkin(12) + qJ(6);
t119 = sin(t123);
t124 = pkin(11) + qJ(4);
t122 = cos(t124);
t139 = t119 * t122;
t121 = cos(t123);
t138 = t121 * t122;
t130 = cos(qJ(2));
t137 = t122 * t130;
t125 = sin(pkin(10));
t126 = sin(pkin(6));
t136 = t125 * t126;
t127 = cos(pkin(10));
t135 = t126 * t127;
t129 = sin(qJ(2));
t134 = t126 * t129;
t133 = t126 * t130;
t128 = cos(pkin(6));
t132 = t128 * t129;
t131 = t128 * t130;
t120 = sin(t124);
t117 = -t125 * t132 + t127 * t130;
t116 = t125 * t131 + t127 * t129;
t115 = t125 * t130 + t127 * t132;
t114 = t125 * t129 - t127 * t131;
t113 = t128 * t120 + t122 * t134;
t112 = -t120 * t134 + t128 * t122;
t111 = t117 * t122 + t120 * t136;
t110 = -t117 * t120 + t122 * t136;
t109 = t115 * t122 - t120 * t135;
t108 = -t115 * t120 - t122 * t135;
t1 = [0, -t116 * t138 + t117 * t119, 0, t110 * t121, 0, -t111 * t119 + t116 * t121; 0, -t114 * t138 + t115 * t119, 0, t108 * t121, 0, -t109 * t119 + t114 * t121; 0 (t119 * t129 + t121 * t137) * t126, 0, t112 * t121, 0, -t113 * t119 - t121 * t133; 0, t116 * t139 + t117 * t121, 0, -t110 * t119, 0, -t111 * t121 - t116 * t119; 0, t114 * t139 + t115 * t121, 0, -t108 * t119, 0, -t109 * t121 - t114 * t119; 0 (-t119 * t137 + t121 * t129) * t126, 0, -t112 * t119, 0, -t113 * t121 + t119 * t133; 0, -t116 * t120, 0, t111, 0, 0; 0, -t114 * t120, 0, t109, 0, 0; 0, t120 * t133, 0, t113, 0, 0;];
JR_rot  = t1;
