% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRP4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:59
% EndTime: 2019-02-26 20:16:59
% DurationCPUTime: 0.07s
% Computational Cost: add. (116->28), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->36)
t133 = qJ(4) + qJ(5);
t131 = sin(t133);
t140 = cos(qJ(3));
t149 = t131 * t140;
t132 = cos(t133);
t148 = t132 * t140;
t135 = sin(pkin(6));
t138 = sin(qJ(3));
t147 = t135 * t138;
t146 = t135 * t140;
t141 = cos(qJ(2));
t145 = t135 * t141;
t137 = cos(pkin(6));
t139 = sin(qJ(2));
t144 = t137 * t139;
t143 = t137 * t141;
t142 = t140 * t141;
t136 = cos(pkin(11));
t134 = sin(pkin(11));
t129 = t137 * t138 + t139 * t146;
t128 = t137 * t140 - t139 * t147;
t127 = -t134 * t144 + t136 * t141;
t126 = t134 * t143 + t136 * t139;
t125 = t134 * t141 + t136 * t144;
t124 = t134 * t139 - t136 * t143;
t123 = t127 * t140 + t134 * t147;
t122 = -t127 * t138 + t134 * t146;
t121 = t125 * t140 - t136 * t147;
t120 = -t125 * t138 - t136 * t146;
t119 = -t129 * t132 + t131 * t145;
t118 = -t129 * t131 - t132 * t145;
t117 = -t123 * t132 - t126 * t131;
t116 = -t123 * t131 + t126 * t132;
t115 = -t121 * t132 - t124 * t131;
t114 = -t121 * t131 + t124 * t132;
t1 = [0, -t126 * t148 + t127 * t131, t122 * t132, t116, t116, 0; 0, -t124 * t148 + t125 * t131, t120 * t132, t114, t114, 0; 0 (t131 * t139 + t132 * t142) * t135, t128 * t132, t118, t118, 0; 0, t126 * t149 + t127 * t132, -t122 * t131, t117, t117, 0; 0, t124 * t149 + t125 * t132, -t120 * t131, t115, t115, 0; 0 (-t131 * t142 + t132 * t139) * t135, -t128 * t131, t119, t119, 0; 0, -t126 * t138, t123, 0, 0, 0; 0, -t124 * t138, t121, 0, 0, 0; 0, t138 * t145, t129, 0, 0, 0;];
JR_rot  = t1;
