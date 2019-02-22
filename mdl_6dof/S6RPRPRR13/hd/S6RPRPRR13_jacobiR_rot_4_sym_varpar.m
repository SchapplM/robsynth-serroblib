% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:38
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR13_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobiR_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:37:55
% EndTime: 2019-02-22 10:37:55
% DurationCPUTime: 0.06s
% Computational Cost: add. (40->18), mult. (122->36), div. (0->0), fcn. (174->10), ass. (0->29)
t124 = cos(pkin(6));
t119 = sin(pkin(12));
t128 = cos(qJ(1));
t133 = t128 * t119;
t122 = cos(pkin(12));
t126 = sin(qJ(1));
t134 = t126 * t122;
t115 = t124 * t133 + t134;
t125 = sin(qJ(3));
t127 = cos(qJ(3));
t132 = t128 * t122;
t135 = t126 * t119;
t114 = -t124 * t132 + t135;
t120 = sin(pkin(7));
t123 = cos(pkin(7));
t121 = sin(pkin(6));
t137 = t121 * t128;
t131 = t114 * t123 + t120 * t137;
t141 = t115 * t125 + t131 * t127;
t139 = t120 * t124;
t138 = t121 * t126;
t136 = t122 * t123;
t116 = -t124 * t134 - t133;
t130 = t116 * t123 + t120 * t138;
t129 = t115 * t127 - t131 * t125;
t117 = -t124 * t135 + t132;
t113 = t117 * t127 + t130 * t125;
t112 = t117 * t125 - t130 * t127;
t1 = [-t114 * t120 + t123 * t137, 0, 0, 0, 0, 0; -t116 * t120 + t123 * t138, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; t129, 0, t112, 0, 0, 0; -t113, 0, t141, 0, 0, 0; 0, 0, -t127 * t139 + (t119 * t125 - t127 * t136) * t121, 0, 0, 0; -t141, 0, t113, 0, 0, 0; t112, 0, t129, 0, 0, 0; 0, 0, t125 * t139 + (t119 * t127 + t125 * t136) * t121, 0, 0, 0;];
JR_rot  = t1;
