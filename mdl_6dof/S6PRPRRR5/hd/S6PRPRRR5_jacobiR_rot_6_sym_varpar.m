% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:15
% EndTime: 2019-02-26 19:56:15
% DurationCPUTime: 0.07s
% Computational Cost: add. (132->32), mult. (214->65), div. (0->0), fcn. (313->10), ass. (0->37)
t140 = sin(pkin(11));
t142 = cos(pkin(11));
t145 = sin(qJ(2));
t143 = cos(pkin(6));
t147 = cos(qJ(2));
t150 = t143 * t147;
t134 = t140 * t150 + t142 * t145;
t139 = qJ(4) + qJ(5);
t137 = sin(t139);
t138 = cos(t139);
t141 = sin(pkin(6));
t154 = t140 * t141;
t125 = t134 * t138 - t137 * t154;
t144 = sin(qJ(6));
t159 = t125 * t144;
t132 = t140 * t145 - t142 * t150;
t153 = t141 * t142;
t127 = t132 * t138 + t137 * t153;
t158 = t127 * t144;
t152 = t141 * t147;
t130 = -t143 * t137 - t138 * t152;
t157 = t130 * t144;
t156 = t137 * t144;
t146 = cos(qJ(6));
t155 = t137 * t146;
t151 = t143 * t145;
t149 = t144 * t145;
t148 = t145 * t146;
t135 = -t140 * t151 + t142 * t147;
t133 = t140 * t147 + t142 * t151;
t131 = -t137 * t152 + t143 * t138;
t129 = t130 * t146;
t128 = t132 * t137 - t138 * t153;
t126 = t134 * t137 + t138 * t154;
t124 = t127 * t146;
t123 = t125 * t146;
t1 = [0, -t134 * t144 + t135 * t155, 0, t123, t123, -t126 * t144 + t135 * t146; 0, -t132 * t144 + t133 * t155, 0, t124, t124, -t128 * t144 + t133 * t146; 0 (t137 * t148 + t144 * t147) * t141, 0, t129, t129, -t131 * t144 + t141 * t148; 0, -t134 * t146 - t135 * t156, 0, -t159, -t159, -t126 * t146 - t135 * t144; 0, -t132 * t146 - t133 * t156, 0, -t158, -t158, -t128 * t146 - t133 * t144; 0 (-t137 * t149 + t146 * t147) * t141, 0, -t157, -t157, -t131 * t146 - t141 * t149; 0, -t135 * t138, 0, t126, t126, 0; 0, -t133 * t138, 0, t128, t128, 0; 0, -t141 * t145 * t138, 0, t131, t131, 0;];
JR_rot  = t1;
