% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP10_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_jacobiR_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:00
% EndTime: 2019-02-26 22:45:01
% DurationCPUTime: 0.10s
% Computational Cost: add. (74->30), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->36)
t131 = cos(pkin(6));
t134 = sin(qJ(2));
t139 = cos(qJ(1));
t141 = t139 * t134;
t135 = sin(qJ(1));
t138 = cos(qJ(2));
t144 = t135 * t138;
t125 = t131 * t141 + t144;
t133 = sin(qJ(3));
t137 = cos(qJ(3));
t130 = sin(pkin(6));
t147 = t130 * t139;
t119 = -t125 * t137 + t133 * t147;
t140 = t139 * t138;
t145 = t135 * t134;
t124 = -t131 * t140 + t145;
t132 = sin(qJ(4));
t136 = cos(qJ(4));
t154 = t119 * t132 + t124 * t136;
t153 = t119 * t136 - t124 * t132;
t150 = t130 * t133;
t149 = t130 * t137;
t148 = t130 * t138;
t146 = t132 * t137;
t143 = t136 * t137;
t142 = t137 * t138;
t117 = -t125 * t133 - t137 * t147;
t127 = -t131 * t145 + t140;
t126 = t131 * t144 + t141;
t123 = t131 * t133 + t134 * t149;
t122 = t131 * t137 - t134 * t150;
t121 = t127 * t137 + t135 * t150;
t120 = t127 * t133 - t135 * t149;
t116 = t121 * t136 + t126 * t132;
t115 = -t121 * t132 + t126 * t136;
t1 = [t153, -t126 * t143 + t127 * t132, -t120 * t136, t115, 0, 0; t116, -t124 * t143 + t125 * t132, t117 * t136, t154, 0, 0; 0 (t132 * t134 + t136 * t142) * t130, t122 * t136, -t123 * t132 - t136 * t148, 0, 0; -t154, t126 * t146 + t127 * t136, t120 * t132, -t116, 0, 0; t115, t124 * t146 + t125 * t136, -t117 * t132, t153, 0, 0; 0 (-t132 * t142 + t134 * t136) * t130, -t122 * t132, -t123 * t136 + t132 * t148, 0, 0; t117, -t126 * t133, t121, 0, 0, 0; t120, -t124 * t133, -t119, 0, 0, 0; 0, t133 * t148, t123, 0, 0, 0;];
JR_rot  = t1;
