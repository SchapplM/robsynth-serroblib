% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP11_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:59
% EndTime: 2019-02-26 22:14:59
% DurationCPUTime: 0.12s
% Computational Cost: add. (71->31), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->36)
t126 = cos(pkin(6));
t129 = sin(qJ(2));
t134 = cos(qJ(1));
t137 = t134 * t129;
t130 = sin(qJ(1));
t133 = cos(qJ(2));
t139 = t130 * t133;
t121 = t126 * t137 + t139;
t128 = sin(qJ(3));
t132 = cos(qJ(3));
t125 = sin(pkin(6));
t144 = t125 * t134;
t113 = t121 * t128 + t132 * t144;
t136 = t134 * t133;
t140 = t130 * t129;
t120 = -t126 * t136 + t140;
t127 = sin(qJ(5));
t131 = cos(qJ(5));
t150 = -t113 * t127 - t120 * t131;
t149 = t113 * t131 - t120 * t127;
t146 = t125 * t128;
t145 = t125 * t132;
t143 = t127 * t128;
t142 = t127 * t133;
t141 = t128 * t131;
t138 = t131 * t133;
t135 = -t121 * t132 + t128 * t144;
t123 = -t126 * t140 + t136;
t122 = t126 * t139 + t137;
t119 = t126 * t128 + t129 * t145;
t118 = -t126 * t132 + t129 * t146;
t117 = t123 * t132 + t130 * t146;
t116 = t123 * t128 - t130 * t145;
t112 = t116 * t127 + t122 * t131;
t111 = t116 * t131 - t122 * t127;
t1 = [t150, -t122 * t143 + t123 * t131, t117 * t127, 0, t111, 0; t112, -t120 * t143 + t121 * t131, -t135 * t127, 0, t149, 0; 0 (t128 * t142 + t129 * t131) * t125, t119 * t127, 0, t118 * t131 + t125 * t142, 0; -t149, -t122 * t141 - t123 * t127, t117 * t131, 0, -t112, 0; t111, -t120 * t141 - t121 * t127, -t135 * t131, 0, t150, 0; 0 (-t127 * t129 + t128 * t138) * t125, t119 * t131, 0, -t118 * t127 + t125 * t138, 0; t135, -t122 * t132, -t116, 0, 0, 0; t117, -t120 * t132, -t113, 0, 0, 0; 0, t133 * t145, -t118, 0, 0, 0;];
JR_rot  = t1;
