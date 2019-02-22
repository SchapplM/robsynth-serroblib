% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:53
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPP2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:52:58
% EndTime: 2019-02-22 09:52:58
% DurationCPUTime: 0.10s
% Computational Cost: add. (54->26), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->29)
t126 = sin(pkin(6));
t130 = sin(qJ(3));
t142 = t126 * t130;
t133 = cos(qJ(3));
t141 = t126 * t133;
t134 = cos(qJ(2));
t140 = t126 * t134;
t128 = cos(pkin(6));
t131 = sin(qJ(2));
t139 = t128 * t131;
t138 = t128 * t134;
t129 = sin(qJ(4));
t137 = t129 * t133;
t132 = cos(qJ(4));
t136 = t132 * t133;
t135 = t133 * t134;
t127 = cos(pkin(10));
t125 = sin(pkin(10));
t123 = t128 * t130 + t131 * t141;
t122 = t128 * t133 - t131 * t142;
t121 = -t125 * t139 + t127 * t134;
t120 = t125 * t138 + t127 * t131;
t119 = t125 * t134 + t127 * t139;
t118 = t125 * t131 - t127 * t138;
t117 = t121 * t133 + t125 * t142;
t116 = -t121 * t130 + t125 * t141;
t115 = t119 * t133 - t127 * t142;
t114 = -t119 * t130 - t127 * t141;
t1 = [0, -t120 * t136 + t121 * t129, t116 * t132, -t117 * t129 + t120 * t132, 0, 0; 0, -t118 * t136 + t119 * t129, t114 * t132, -t115 * t129 + t118 * t132, 0, 0; 0 (t129 * t131 + t132 * t135) * t126, t122 * t132, -t123 * t129 - t132 * t140, 0, 0; 0, -t120 * t137 - t121 * t132, t116 * t129, t117 * t132 + t120 * t129, 0, 0; 0, -t118 * t137 - t119 * t132, t114 * t129, t115 * t132 + t118 * t129, 0, 0; 0 (t129 * t135 - t131 * t132) * t126, t122 * t129, t123 * t132 - t129 * t140, 0, 0; 0, t120 * t130, -t117, 0, 0, 0; 0, t118 * t130, -t115, 0, 0, 0; 0, -t130 * t140, -t123, 0, 0, 0;];
JR_rot  = t1;
