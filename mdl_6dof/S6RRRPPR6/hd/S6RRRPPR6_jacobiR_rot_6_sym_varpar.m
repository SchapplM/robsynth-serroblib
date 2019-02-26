% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:06:23
% EndTime: 2019-02-26 22:06:23
% DurationCPUTime: 0.17s
% Computational Cost: add. (122->32), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
t138 = cos(pkin(6));
t140 = sin(qJ(2));
t144 = cos(qJ(1));
t147 = t144 * t140;
t141 = sin(qJ(1));
t143 = cos(qJ(2));
t149 = t141 * t143;
t130 = t138 * t147 + t149;
t136 = qJ(3) + pkin(11);
t134 = sin(t136);
t135 = cos(t136);
t137 = sin(pkin(6));
t152 = t137 * t144;
t122 = t130 * t134 + t135 * t152;
t146 = t144 * t143;
t150 = t141 * t140;
t129 = -t138 * t146 + t150;
t139 = sin(qJ(6));
t142 = cos(qJ(6));
t160 = -t122 * t139 - t129 * t142;
t159 = t122 * t142 - t129 * t139;
t156 = t134 * t139;
t155 = t134 * t142;
t154 = t137 * t140;
t153 = t137 * t141;
t151 = t139 * t143;
t148 = t142 * t143;
t145 = -t130 * t135 + t134 * t152;
t132 = -t138 * t150 + t146;
t131 = t138 * t149 + t147;
t128 = t138 * t134 + t135 * t154;
t127 = t134 * t154 - t138 * t135;
t126 = t132 * t135 + t134 * t153;
t125 = t132 * t134 - t135 * t153;
t121 = t125 * t139 + t131 * t142;
t120 = t125 * t142 - t131 * t139;
t1 = [t160, -t131 * t156 + t132 * t142, t126 * t139, 0, 0, t120; t121, -t129 * t156 + t130 * t142, -t145 * t139, 0, 0, t159; 0 (t134 * t151 + t140 * t142) * t137, t128 * t139, 0, 0, t127 * t142 + t137 * t151; -t159, -t131 * t155 - t132 * t139, t126 * t142, 0, 0, -t121; t120, -t129 * t155 - t130 * t139, -t145 * t142, 0, 0, t160; 0 (t134 * t148 - t139 * t140) * t137, t128 * t142, 0, 0, -t127 * t139 + t137 * t148; t145, -t131 * t135, -t125, 0, 0, 0; t126, -t129 * t135, -t122, 0, 0, 0; 0, t137 * t143 * t135, -t127, 0, 0, 0;];
JR_rot  = t1;
