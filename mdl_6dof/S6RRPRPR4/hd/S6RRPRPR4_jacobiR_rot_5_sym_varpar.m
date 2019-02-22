% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:24
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:24:31
% EndTime: 2019-02-22 11:24:31
% DurationCPUTime: 0.06s
% Computational Cost: add. (92->19), mult. (182->36), div. (0->0), fcn. (266->10), ass. (0->28)
t127 = sin(pkin(6));
t131 = sin(qJ(1));
t138 = t127 * t131;
t133 = cos(qJ(1));
t137 = t127 * t133;
t129 = cos(pkin(6));
t126 = sin(pkin(11));
t128 = cos(pkin(11));
t130 = sin(qJ(2));
t132 = cos(qJ(2));
t135 = t132 * t126 + t130 * t128;
t118 = t135 * t129;
t119 = t130 * t126 - t132 * t128;
t112 = t133 * t118 - t131 * t119;
t125 = qJ(4) + pkin(12);
t123 = sin(t125);
t124 = cos(t125);
t136 = -t112 * t124 + t123 * t137;
t114 = -t131 * t118 - t133 * t119;
t134 = t112 * t123 + t124 * t137;
t117 = t119 * t129;
t116 = t135 * t127;
t115 = t119 * t127;
t113 = t131 * t117 - t133 * t135;
t111 = -t133 * t117 - t131 * t135;
t110 = t114 * t124 + t123 * t138;
t109 = -t114 * t123 + t124 * t138;
t1 = [t136, t113 * t124, 0, t109, 0, 0; t110, t111 * t124, 0, -t134, 0, 0; 0, -t115 * t124, 0, -t116 * t123 + t129 * t124, 0, 0; t134, -t113 * t123, 0, -t110, 0, 0; t109, -t111 * t123, 0, t136, 0, 0; 0, t115 * t123, 0, -t116 * t124 - t129 * t123, 0, 0; t111, t114, 0, 0, 0, 0; -t113, t112, 0, 0, 0, 0; 0, t116, 0, 0, 0, 0;];
JR_rot  = t1;
