% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR5
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
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:17
% EndTime: 2019-02-26 22:18:17
% DurationCPUTime: 0.04s
% Computational Cost: add. (95->16), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
t131 = qJ(5) + qJ(6);
t129 = cos(t131);
t132 = qJ(2) + qJ(3);
t130 = cos(t132);
t141 = t130 * t129;
t128 = sin(t132);
t133 = sin(qJ(1));
t140 = t133 * t128;
t139 = t133 * t129;
t138 = t133 * t130;
t134 = cos(qJ(1));
t137 = t134 * t128;
t136 = t134 * t129;
t135 = t134 * t130;
t127 = sin(t131);
t126 = t130 * t127;
t125 = t128 * t129;
t124 = t128 * t127;
t123 = t129 * t135;
t122 = t127 * t135;
t121 = t129 * t138;
t120 = t127 * t138;
t119 = -t127 * t140 + t136;
t118 = t134 * t127 + t128 * t139;
t117 = t127 * t137 + t139;
t116 = -t133 * t127 + t128 * t136;
t1 = [t119, t122, t122, 0, t116, t116; t117, t120, t120, 0, t118, t118; 0, t124, t124, 0, -t141, -t141; -t118, t123, t123, 0, -t117, -t117; t116, t121, t121, 0, t119, t119; 0, t125, t125, 0, t126, t126; -t138, -t137, -t137, 0, 0, 0; t135, -t140, -t140, 0, 0, 0; 0, t130, t130, 0, 0, 0;];
JR_rot  = t1;
