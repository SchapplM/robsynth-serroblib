% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR14_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:26
% EndTime: 2019-02-26 21:45:26
% DurationCPUTime: 0.12s
% Computational Cost: add. (71->28), mult. (219->61), div. (0->0), fcn. (320->10), ass. (0->35)
t130 = cos(pkin(6));
t133 = sin(qJ(2));
t138 = cos(qJ(1));
t142 = t138 * t133;
t134 = sin(qJ(1));
t137 = cos(qJ(2));
t144 = t134 * t137;
t125 = t130 * t142 + t144;
t131 = sin(qJ(6));
t135 = cos(qJ(6));
t141 = t138 * t137;
t145 = t134 * t133;
t124 = -t130 * t141 + t145;
t132 = sin(qJ(4));
t136 = cos(qJ(4));
t129 = sin(pkin(6));
t147 = t129 * t138;
t139 = t124 * t136 + t132 * t147;
t154 = -t125 * t135 + t131 * t139;
t153 = t125 * t131 + t135 * t139;
t150 = t129 * t133;
t149 = t129 * t134;
t148 = t129 * t137;
t146 = t131 * t136;
t143 = t135 * t136;
t140 = -t124 * t132 + t136 * t147;
t127 = -t130 * t145 + t141;
t126 = t130 * t144 + t142;
t123 = t130 * t136 - t132 * t148;
t122 = t130 * t132 + t136 * t148;
t118 = t126 * t132 + t136 * t149;
t117 = -t126 * t136 + t132 * t149;
t116 = t117 * t131 + t127 * t135;
t115 = t117 * t135 - t127 * t131;
t1 = [t154, -t126 * t135 - t127 * t146, 0, t118 * t131, 0, t115; t116, -t124 * t135 - t125 * t146, 0, -t140 * t131, 0, -t153; 0 (-t133 * t146 + t135 * t137) * t129, 0, t123 * t131, 0, t122 * t135 - t131 * t150; t153, t126 * t131 - t127 * t143, 0, t118 * t135, 0, -t116; t115, t124 * t131 - t125 * t143, 0, -t140 * t135, 0, t154; 0 (-t131 * t137 - t133 * t143) * t129, 0, t123 * t135, 0, -t122 * t131 - t135 * t150; t140, t127 * t132, 0, -t117, 0, 0; t118, t125 * t132, 0, t139, 0, 0; 0, t132 * t150, 0, -t122, 0, 0;];
JR_rot  = t1;
