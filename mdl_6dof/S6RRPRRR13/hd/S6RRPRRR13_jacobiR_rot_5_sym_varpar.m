% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR13_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:00:59
% EndTime: 2019-02-26 22:01:00
% DurationCPUTime: 0.13s
% Computational Cost: add. (77->31), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->36)
t131 = cos(pkin(6));
t138 = cos(qJ(2));
t139 = cos(qJ(1));
t140 = t139 * t138;
t134 = sin(qJ(2));
t135 = sin(qJ(1));
t143 = t135 * t134;
t124 = -t131 * t140 + t143;
t133 = sin(qJ(4));
t137 = cos(qJ(4));
t130 = sin(pkin(6));
t148 = t130 * t139;
t119 = -t124 * t133 + t137 * t148;
t141 = t139 * t134;
t142 = t135 * t138;
t125 = t131 * t141 + t142;
t132 = sin(qJ(5));
t136 = cos(qJ(5));
t154 = t119 * t132 + t125 * t136;
t153 = t119 * t136 - t125 * t132;
t150 = t130 * t137;
t149 = t130 * t138;
t147 = t132 * t133;
t146 = t132 * t134;
t145 = t133 * t136;
t144 = t134 * t136;
t118 = t124 * t137 + t133 * t148;
t127 = -t131 * t143 + t140;
t126 = t131 * t142 + t141;
t123 = t131 * t137 - t133 * t149;
t122 = -t131 * t133 - t137 * t149;
t117 = t126 * t133 + t135 * t150;
t116 = t135 * t130 * t133 - t126 * t137;
t115 = t117 * t136 + t127 * t132;
t114 = -t117 * t132 + t127 * t136;
t1 = [t153, -t126 * t132 + t127 * t145, 0, -t116 * t136, t114, 0; t115, -t124 * t132 + t125 * t145, 0, t118 * t136, t154, 0; 0 (t132 * t138 + t133 * t144) * t130, 0, t122 * t136, -t123 * t132 + t130 * t144, 0; -t154, -t126 * t136 - t127 * t147, 0, t116 * t132, -t115, 0; t114, -t124 * t136 - t125 * t147, 0, -t118 * t132, t153, 0; 0 (-t133 * t146 + t136 * t138) * t130, 0, -t122 * t132, -t123 * t136 - t130 * t146, 0; t118, -t127 * t137, 0, t117, 0, 0; t116, -t125 * t137, 0, -t119, 0, 0; 0, -t134 * t150, 0, t123, 0, 0;];
JR_rot  = t1;
