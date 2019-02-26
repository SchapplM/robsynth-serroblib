% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:37
% EndTime: 2019-02-26 20:07:37
% DurationCPUTime: 0.07s
% Computational Cost: add. (113->28), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->36)
t135 = qJ(5) + qJ(6);
t133 = sin(t135);
t140 = sin(qJ(3));
t151 = t133 * t140;
t134 = cos(t135);
t150 = t134 * t140;
t137 = sin(pkin(6));
t149 = t137 * t140;
t142 = cos(qJ(3));
t148 = t137 * t142;
t143 = cos(qJ(2));
t147 = t137 * t143;
t139 = cos(pkin(6));
t141 = sin(qJ(2));
t146 = t139 * t141;
t145 = t139 * t143;
t144 = t140 * t143;
t138 = cos(pkin(11));
t136 = sin(pkin(11));
t131 = t139 * t140 + t141 * t148;
t130 = -t139 * t142 + t141 * t149;
t129 = -t136 * t146 + t138 * t143;
t128 = t136 * t145 + t138 * t141;
t127 = t136 * t143 + t138 * t146;
t126 = t136 * t141 - t138 * t145;
t125 = t129 * t142 + t136 * t149;
t124 = t129 * t140 - t136 * t148;
t123 = t127 * t142 - t138 * t149;
t122 = t127 * t140 + t138 * t148;
t121 = -t130 * t133 + t134 * t147;
t120 = t130 * t134 + t133 * t147;
t119 = -t124 * t133 - t128 * t134;
t118 = t124 * t134 - t128 * t133;
t117 = -t122 * t133 - t126 * t134;
t116 = t122 * t134 - t126 * t133;
t1 = [0, -t128 * t151 + t129 * t134, t125 * t133, 0, t118, t118; 0, -t126 * t151 + t127 * t134, t123 * t133, 0, t116, t116; 0 (t133 * t144 + t134 * t141) * t137, t131 * t133, 0, t120, t120; 0, -t128 * t150 - t129 * t133, t125 * t134, 0, t119, t119; 0, -t126 * t150 - t127 * t133, t123 * t134, 0, t117, t117; 0 (-t133 * t141 + t134 * t144) * t137, t131 * t134, 0, t121, t121; 0, -t128 * t142, -t124, 0, 0, 0; 0, -t126 * t142, -t122, 0, 0, 0; 0, t142 * t147, -t130, 0, 0, 0;];
JR_rot  = t1;
