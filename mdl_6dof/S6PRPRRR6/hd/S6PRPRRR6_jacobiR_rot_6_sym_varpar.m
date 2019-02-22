% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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
% Datum: 2019-02-22 09:40
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:40:45
% EndTime: 2019-02-22 09:40:45
% DurationCPUTime: 0.11s
% Computational Cost: add. (119->29), mult. (219->64), div. (0->0), fcn. (320->10), ass. (0->37)
t133 = qJ(5) + qJ(6);
t131 = sin(t133);
t138 = sin(qJ(4));
t150 = t131 * t138;
t132 = cos(t133);
t149 = t132 * t138;
t135 = sin(pkin(6));
t148 = t135 * t138;
t139 = sin(qJ(2));
t147 = t135 * t139;
t140 = cos(qJ(4));
t146 = t135 * t140;
t141 = cos(qJ(2));
t145 = t135 * t141;
t137 = cos(pkin(6));
t144 = t137 * t139;
t143 = t137 * t141;
t142 = t138 * t139;
t136 = cos(pkin(11));
t134 = sin(pkin(11));
t129 = t137 * t140 - t138 * t145;
t128 = -t137 * t138 - t140 * t145;
t127 = -t134 * t144 + t136 * t141;
t126 = t134 * t143 + t136 * t139;
t125 = t134 * t141 + t136 * t144;
t124 = t134 * t139 - t136 * t143;
t123 = t124 * t138 - t136 * t146;
t122 = t124 * t140 + t136 * t148;
t121 = t126 * t138 + t134 * t146;
t120 = t126 * t140 - t134 * t148;
t119 = -t129 * t132 - t131 * t147;
t118 = -t129 * t131 + t132 * t147;
t117 = -t123 * t132 - t125 * t131;
t116 = -t123 * t131 + t125 * t132;
t115 = -t121 * t132 - t127 * t131;
t114 = -t121 * t131 + t127 * t132;
t1 = [0, -t126 * t131 + t127 * t149, 0, t120 * t132, t114, t114; 0, -t124 * t131 + t125 * t149, 0, t122 * t132, t116, t116; 0 (t131 * t141 + t132 * t142) * t135, 0, t128 * t132, t118, t118; 0, -t126 * t132 - t127 * t150, 0, -t120 * t131, t115, t115; 0, -t124 * t132 - t125 * t150, 0, -t122 * t131, t117, t117; 0 (-t131 * t142 + t132 * t141) * t135, 0, -t128 * t131, t119, t119; 0, -t127 * t140, 0, t121, 0, 0; 0, -t125 * t140, 0, t123, 0, 0; 0, -t139 * t146, 0, t129, 0, 0;];
JR_rot  = t1;
