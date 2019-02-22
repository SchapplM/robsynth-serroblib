% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:22
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR10_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:22:47
% EndTime: 2019-02-22 12:22:47
% DurationCPUTime: 0.05s
% Computational Cost: add. (77->16), mult. (117->30), div. (0->0), fcn. (180->8), ass. (0->28)
t130 = sin(pkin(6));
t132 = sin(qJ(2));
t143 = t130 * t132;
t133 = sin(qJ(1));
t142 = t130 * t133;
t134 = cos(qJ(2));
t141 = t130 * t134;
t135 = cos(qJ(1));
t140 = t130 * t135;
t139 = t133 * t132;
t138 = t133 * t134;
t137 = t135 * t132;
t136 = t135 * t134;
t131 = cos(pkin(6));
t124 = t131 * t137 + t138;
t129 = qJ(3) + qJ(4);
t127 = sin(t129);
t128 = cos(t129);
t117 = t124 * t127 + t128 * t140;
t118 = t124 * t128 - t127 * t140;
t126 = -t131 * t139 + t136;
t125 = t131 * t138 + t137;
t123 = t131 * t136 - t139;
t122 = t131 * t127 + t128 * t143;
t121 = t127 * t143 - t131 * t128;
t120 = t126 * t128 + t127 * t142;
t119 = t126 * t127 - t128 * t142;
t1 = [t123, t126, 0, 0, 0, 0; t125, t124, 0, 0, 0, 0; 0, t143, 0, 0, 0, 0; t118, t125 * t128, t119, t119, 0, 0; -t120, -t123 * t128, t117, t117, 0, 0; 0, -t128 * t141, t121, t121, 0, 0; -t117, -t125 * t127, t120, t120, 0, 0; t119, t123 * t127, t118, t118, 0, 0; 0, t127 * t141, t122, t122, 0, 0;];
JR_rot  = t1;
