% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:55:23
% EndTime: 2019-02-26 19:55:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (167->29), mult. (219->64), div. (0->0), fcn. (320->10), ass. (0->38)
t142 = pkin(12) + qJ(4);
t139 = cos(t142);
t143 = qJ(5) + qJ(6);
t140 = sin(t143);
t158 = t139 * t140;
t141 = cos(t143);
t157 = t139 * t141;
t149 = cos(qJ(2));
t156 = t139 * t149;
t144 = sin(pkin(11));
t145 = sin(pkin(6));
t155 = t144 * t145;
t146 = cos(pkin(11));
t154 = t145 * t146;
t148 = sin(qJ(2));
t153 = t145 * t148;
t152 = t145 * t149;
t147 = cos(pkin(6));
t151 = t147 * t148;
t150 = t147 * t149;
t138 = sin(t142);
t136 = -t144 * t151 + t146 * t149;
t135 = t144 * t150 + t146 * t148;
t134 = t144 * t149 + t146 * t151;
t133 = t144 * t148 - t146 * t150;
t132 = t147 * t138 + t139 * t153;
t131 = -t138 * t153 + t147 * t139;
t130 = t136 * t139 + t138 * t155;
t129 = -t136 * t138 + t139 * t155;
t128 = t134 * t139 - t138 * t154;
t127 = -t134 * t138 - t139 * t154;
t126 = -t132 * t141 + t140 * t152;
t125 = -t132 * t140 - t141 * t152;
t124 = -t130 * t141 - t135 * t140;
t123 = -t130 * t140 + t135 * t141;
t122 = -t128 * t141 - t133 * t140;
t121 = -t128 * t140 + t133 * t141;
t1 = [0, -t135 * t157 + t136 * t140, 0, t129 * t141, t123, t123; 0, -t133 * t157 + t134 * t140, 0, t127 * t141, t121, t121; 0 (t140 * t148 + t141 * t156) * t145, 0, t131 * t141, t125, t125; 0, t135 * t158 + t136 * t141, 0, -t129 * t140, t124, t124; 0, t133 * t158 + t134 * t141, 0, -t127 * t140, t122, t122; 0 (-t140 * t156 + t141 * t148) * t145, 0, -t131 * t140, t126, t126; 0, -t135 * t138, 0, t130, 0, 0; 0, -t133 * t138, 0, t128, 0, 0; 0, t138 * t152, 0, t132, 0, 0;];
JR_rot  = t1;
