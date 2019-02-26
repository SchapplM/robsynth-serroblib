% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:37
% EndTime: 2019-02-26 21:40:37
% DurationCPUTime: 0.06s
% Computational Cost: add. (66->18), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->27)
t127 = sin(pkin(6));
t132 = sin(qJ(1));
t140 = t127 * t132;
t135 = cos(qJ(1));
t139 = t127 * t135;
t129 = cos(pkin(6));
t126 = sin(pkin(11));
t128 = cos(pkin(11));
t131 = sin(qJ(2));
t134 = cos(qJ(2));
t138 = t134 * t126 + t131 * t128;
t122 = t138 * t129;
t123 = t131 * t126 - t134 * t128;
t116 = t135 * t122 - t132 * t123;
t118 = -t132 * t122 - t135 * t123;
t130 = sin(qJ(4));
t133 = cos(qJ(4));
t137 = t116 * t130 + t133 * t139;
t136 = t116 * t133 - t130 * t139;
t121 = t123 * t129;
t120 = t138 * t127;
t119 = t123 * t127;
t117 = t132 * t121 - t135 * t138;
t115 = -t135 * t121 - t132 * t138;
t114 = t118 * t133 + t130 * t140;
t113 = t118 * t130 - t133 * t140;
t1 = [t115, t118, 0, 0, 0, 0; -t117, t116, 0, 0, 0, 0; 0, t120, 0, 0, 0, 0; t136, -t117 * t133, 0, t113, 0, 0; -t114, -t115 * t133, 0, t137, 0, 0; 0, t119 * t133, 0, t120 * t130 - t129 * t133, 0, 0; -t137, t117 * t130, 0, t114, 0, 0; t113, t115 * t130, 0, t136, 0, 0; 0, -t119 * t130, 0, t120 * t133 + t129 * t130, 0, 0;];
JR_rot  = t1;
