% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:52
% EndTime: 2019-02-26 20:19:52
% DurationCPUTime: 0.08s
% Computational Cost: add. (202->28), mult. (275->63), div. (0->0), fcn. (402->10), ass. (0->36)
t138 = qJ(4) + qJ(5) + qJ(6);
t136 = sin(t138);
t145 = cos(qJ(3));
t154 = t136 * t145;
t137 = cos(t138);
t153 = t137 * t145;
t140 = sin(pkin(6));
t143 = sin(qJ(3));
t152 = t140 * t143;
t151 = t140 * t145;
t146 = cos(qJ(2));
t150 = t140 * t146;
t142 = cos(pkin(6));
t144 = sin(qJ(2));
t149 = t142 * t144;
t148 = t142 * t146;
t147 = t145 * t146;
t141 = cos(pkin(12));
t139 = sin(pkin(12));
t134 = t142 * t143 + t144 * t151;
t133 = t142 * t145 - t144 * t152;
t132 = -t139 * t149 + t141 * t146;
t131 = t139 * t148 + t141 * t144;
t130 = t139 * t146 + t141 * t149;
t129 = t139 * t144 - t141 * t148;
t128 = t132 * t145 + t139 * t152;
t127 = -t132 * t143 + t139 * t151;
t126 = t130 * t145 - t141 * t152;
t125 = -t130 * t143 - t141 * t151;
t124 = -t134 * t137 + t136 * t150;
t123 = -t134 * t136 - t137 * t150;
t122 = -t128 * t137 - t131 * t136;
t121 = -t128 * t136 + t131 * t137;
t120 = -t126 * t137 - t129 * t136;
t119 = -t126 * t136 + t129 * t137;
t1 = [0, -t131 * t153 + t132 * t136, t127 * t137, t121, t121, t121; 0, -t129 * t153 + t130 * t136, t125 * t137, t119, t119, t119; 0 (t136 * t144 + t137 * t147) * t140, t133 * t137, t123, t123, t123; 0, t131 * t154 + t132 * t137, -t127 * t136, t122, t122, t122; 0, t129 * t154 + t130 * t137, -t125 * t136, t120, t120, t120; 0 (-t136 * t147 + t137 * t144) * t140, -t133 * t136, t124, t124, t124; 0, -t131 * t143, t128, 0, 0, 0; 0, -t129 * t143, t126, 0, 0, 0; 0, t143 * t150, t134, 0, 0, 0;];
JR_rot  = t1;
