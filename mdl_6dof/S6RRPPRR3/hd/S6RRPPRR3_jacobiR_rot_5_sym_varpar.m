% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:13
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:13:43
% EndTime: 2019-02-22 11:13:43
% DurationCPUTime: 0.06s
% Computational Cost: add. (92->19), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->28)
t126 = sin(pkin(6));
t130 = sin(qJ(1));
t137 = t126 * t130;
t132 = cos(qJ(1));
t136 = t126 * t132;
t128 = cos(pkin(6));
t125 = sin(pkin(11));
t127 = cos(pkin(11));
t129 = sin(qJ(2));
t131 = cos(qJ(2));
t134 = t131 * t125 + t129 * t127;
t117 = t134 * t128;
t118 = t129 * t125 - t131 * t127;
t111 = t132 * t117 - t130 * t118;
t124 = pkin(12) + qJ(5);
t122 = sin(t124);
t123 = cos(t124);
t135 = -t111 * t123 + t122 * t136;
t113 = -t130 * t117 - t132 * t118;
t133 = t111 * t122 + t123 * t136;
t116 = t118 * t128;
t115 = t134 * t126;
t114 = t118 * t126;
t112 = t130 * t116 - t132 * t134;
t110 = -t132 * t116 - t130 * t134;
t109 = t113 * t123 + t122 * t137;
t108 = -t113 * t122 + t123 * t137;
t1 = [t135, t112 * t123, 0, 0, t108, 0; t109, t110 * t123, 0, 0, -t133, 0; 0, -t114 * t123, 0, 0, -t115 * t122 + t128 * t123, 0; t133, -t112 * t122, 0, 0, -t109, 0; t108, -t110 * t122, 0, 0, t135, 0; 0, t114 * t122, 0, 0, -t115 * t123 - t128 * t122, 0; t110, t113, 0, 0, 0, 0; -t112, t111, 0, 0, 0, 0; 0, t115, 0, 0, 0, 0;];
JR_rot  = t1;
