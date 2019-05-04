% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR10V2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiR_rot_5_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:32
% EndTime: 2019-04-11 14:56:32
% DurationCPUTime: 0.10s
% Computational Cost: add. (106->23), mult. (149->41), div. (0->0), fcn. (226->8), ass. (0->38)
t131 = qJ(2) + qJ(3);
t129 = sin(t131);
t132 = sin(qJ(5));
t152 = t129 * t132;
t135 = cos(qJ(5));
t151 = t129 * t135;
t137 = cos(qJ(1));
t150 = t129 * t137;
t136 = cos(qJ(4));
t149 = t132 * t136;
t133 = sin(qJ(4));
t134 = sin(qJ(1));
t148 = t134 * t133;
t147 = t134 * t136;
t146 = t135 * t136;
t145 = t137 * t133;
t144 = t137 * t136;
t143 = t129 * t148;
t142 = t129 * t145;
t130 = cos(t131);
t123 = t130 * t147 - t145;
t141 = -t123 * t132 + t134 * t151;
t140 = -t123 * t135 - t134 * t152;
t139 = -t129 * t146 + t130 * t132;
t138 = t129 * t149 + t130 * t135;
t126 = t130 * t133;
t125 = t130 * t144 + t148;
t124 = t130 * t145 - t147;
t122 = -t130 * t148 - t144;
t121 = t130 * t146 + t152;
t120 = -t130 * t149 + t151;
t119 = t139 * t137;
t118 = t138 * t137;
t117 = t139 * t134;
t116 = t138 * t134;
t115 = t125 * t135 + t132 * t150;
t114 = -t125 * t132 + t135 * t150;
t1 = [t140, t119, t119, -t124 * t135, t114, 0; t115, t117, t117, t122 * t135, t141, 0; 0, t121, t121, -t133 * t151, -t138, 0; -t141, t118, t118, t124 * t132, -t115, 0; t114, t116, t116, -t122 * t132, t140, 0; 0, t120, t120, t133 * t152, t139, 0; t122, -t142, -t142, t125, 0, 0; t124, -t143, -t143, t123, 0, 0; 0, t126, t126, t129 * t136, 0, 0;];
JR_rot  = t1;
