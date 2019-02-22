% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:54
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:54:07
% EndTime: 2019-02-22 11:54:07
% DurationCPUTime: 0.13s
% Computational Cost: add. (74->31), mult. (219->64), div. (0->0), fcn. (320->10), ass. (0->34)
t134 = sin(qJ(2));
t135 = sin(qJ(1));
t138 = cos(qJ(2));
t139 = cos(qJ(1));
t152 = cos(pkin(6));
t141 = t139 * t152;
t126 = t134 * t141 + t135 * t138;
t133 = sin(qJ(3));
t137 = cos(qJ(3));
t131 = sin(pkin(6));
t147 = t131 * t139;
t118 = t126 * t133 + t137 * t147;
t125 = t135 * t134 - t138 * t141;
t132 = sin(qJ(6));
t136 = cos(qJ(6));
t154 = t118 * t132 + t125 * t136;
t153 = -t118 * t136 + t125 * t132;
t149 = t131 * t133;
t148 = t131 * t137;
t146 = t132 * t133;
t145 = t132 * t138;
t144 = t133 * t136;
t143 = t136 * t138;
t142 = t135 * t152;
t140 = -t126 * t137 + t133 * t147;
t128 = -t134 * t142 + t139 * t138;
t127 = -t139 * t134 - t138 * t142;
t124 = t152 * t133 + t134 * t148;
t123 = t134 * t149 - t152 * t137;
t122 = t128 * t137 + t135 * t149;
t121 = t128 * t133 - t135 * t148;
t117 = t121 * t136 + t127 * t132;
t116 = -t121 * t132 + t127 * t136;
t1 = [t153, t127 * t144 - t128 * t132, t122 * t136, 0, 0, t116; t117, -t125 * t144 - t126 * t132, -t140 * t136, 0, 0, -t154; 0 (-t132 * t134 + t133 * t143) * t131, t124 * t136, 0, 0, -t123 * t132 + t131 * t143; t154, -t127 * t146 - t128 * t136, -t122 * t132, 0, 0, -t117; t116, t125 * t146 - t126 * t136, t140 * t132, 0, 0, t153; 0 (-t133 * t145 - t134 * t136) * t131, -t124 * t132, 0, 0, -t123 * t136 - t131 * t145; t140, t127 * t137, -t121, 0, 0, 0; t122, -t125 * t137, -t118, 0, 0, 0; 0, t138 * t148, -t123, 0, 0, 0;];
JR_rot  = t1;
