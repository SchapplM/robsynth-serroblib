% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR9_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiR_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:27
% EndTime: 2019-02-26 22:20:28
% DurationCPUTime: 0.11s
% Computational Cost: add. (100->28), mult. (288->62), div. (0->0), fcn. (407->12), ass. (0->31)
t146 = cos(pkin(6));
t129 = sin(pkin(6));
t134 = sin(qJ(1));
t145 = t129 * t134;
t137 = cos(qJ(1));
t144 = t129 * t137;
t143 = t134 * t146;
t142 = t137 * t146;
t127 = sin(pkin(13));
t130 = cos(pkin(13));
t132 = sin(qJ(3));
t135 = cos(qJ(3));
t141 = t135 * t127 + t132 * t130;
t122 = t132 * t127 - t135 * t130;
t128 = sin(pkin(7));
t114 = t122 * t128;
t131 = cos(pkin(7));
t116 = t122 * t131;
t133 = sin(qJ(2));
t136 = cos(qJ(2));
t118 = t134 * t133 - t136 * t142;
t119 = t133 * t142 + t134 * t136;
t140 = -t114 * t144 - t118 * t116 + t119 * t141;
t115 = t141 * t128;
t117 = t141 * t131;
t139 = t115 * t144 + t118 * t117 + t119 * t122;
t120 = -t137 * t133 - t136 * t143;
t121 = -t133 * t143 + t137 * t136;
t138 = t115 * t145 + t120 * t117 - t121 * t122;
t113 = -t114 * t145 - t120 * t116 - t121 * t141;
t1 = [t139, -t121 * t117 - t120 * t122, t113, 0, 0, 0; t138, -t119 * t117 + t118 * t122, -t140, 0, 0, 0; 0 (-t117 * t133 - t122 * t136) * t129, -t146 * t114 + (-t116 * t136 - t133 * t141) * t129, 0, 0, 0; t140, t121 * t116 - t120 * t141, -t138, 0, 0, 0; t113, t119 * t116 + t118 * t141, t139, 0, 0, 0; 0 (t116 * t133 - t136 * t141) * t129, -t146 * t115 + (-t117 * t136 + t122 * t133) * t129, 0, 0, 0; -t118 * t128 + t131 * t144, t121 * t128, 0, 0, 0, 0; -t120 * t128 + t131 * t145, t119 * t128, 0, 0, 0, 0; 0, t129 * t133 * t128, 0, 0, 0, 0;];
JR_rot  = t1;
