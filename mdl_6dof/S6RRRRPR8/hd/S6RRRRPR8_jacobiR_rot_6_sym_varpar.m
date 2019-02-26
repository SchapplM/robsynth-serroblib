% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:43
% EndTime: 2019-02-26 22:34:43
% DurationCPUTime: 0.13s
% Computational Cost: add. (152->30), mult. (188->32), div. (0->0), fcn. (281->8), ass. (0->24)
t131 = sin(qJ(2));
t129 = qJ(3) + qJ(4);
t127 = sin(t129);
t128 = cos(t129);
t130 = sin(qJ(6));
t133 = cos(qJ(6));
t137 = t127 * t130 + t128 * t133;
t144 = t137 * t131;
t135 = cos(qJ(1));
t132 = sin(qJ(1));
t134 = cos(qJ(2));
t141 = t132 * t134;
t121 = t127 * t141 + t135 * t128;
t122 = -t135 * t127 + t128 * t141;
t143 = t121 * t133 - t122 * t130;
t140 = t135 * t134;
t112 = t121 * t130 + t122 * t133;
t123 = t127 * t140 - t132 * t128;
t124 = t132 * t127 + t128 * t140;
t116 = t123 * t133 - t124 * t130;
t139 = t123 * t130 + t124 * t133;
t138 = t127 * t133 - t128 * t130;
t119 = t138 * t131;
t1 = [-t112, -t135 * t144, -t116, -t116, 0, t116; t139, -t132 * t144, -t143, -t143, 0, t143; 0, t137 * t134, -t119, -t119, 0, t119; -t143, -t135 * t119, t139, t139, 0, -t139; t116, -t132 * t119, t112, t112, 0, -t112; 0, t138 * t134, t144, t144, 0, -t144; t132 * t131, -t140, 0, 0, 0, 0; -t135 * t131, -t141, 0, 0, 0, 0; 0, -t131, 0, 0, 0, 0;];
JR_rot  = t1;
