% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:10
% EndTime: 2019-02-26 22:41:10
% DurationCPUTime: 0.04s
% Computational Cost: add. (99->20), mult. (66->20), div. (0->0), fcn. (114->6), ass. (0->27)
t134 = qJ(4) + qJ(5);
t130 = sin(t134);
t135 = qJ(2) + qJ(3);
t131 = sin(t135);
t146 = t131 * t130;
t136 = sin(qJ(1));
t145 = t136 * t131;
t132 = cos(t134);
t144 = t136 * t132;
t133 = cos(t135);
t128 = t136 * t133;
t137 = cos(qJ(1));
t143 = t137 * t131;
t142 = t137 * t132;
t129 = t137 * t133;
t141 = t130 * t145;
t140 = t131 * t144;
t139 = t130 * t143;
t138 = t131 * t142;
t126 = t133 * t132;
t125 = t133 * t130;
t124 = t131 * t132;
t123 = t132 * t129 + t136 * t130;
t122 = t130 * t129 - t144;
t121 = t132 * t128 - t137 * t130;
t120 = -t130 * t128 - t142;
t1 = [-t121, -t138, -t138, -t122, -t122, 0; t123, -t140, -t140, t120, t120, 0; 0, t126, t126, -t146, -t146, 0; -t145, t129, t129, 0, 0, 0; t143, t128, t128, 0, 0, 0; 0, t131, t131, 0, 0, 0; t120, -t139, -t139, t123, t123, 0; t122, -t141, -t141, t121, t121, 0; 0, t125, t125, t124, t124, 0;];
JR_rot  = t1;
