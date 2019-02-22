% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:50
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:50:15
% EndTime: 2019-02-22 09:50:15
% DurationCPUTime: 0.07s
% Computational Cost: add. (158->28), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->36)
t137 = pkin(12) + qJ(5) + qJ(6);
t135 = sin(t137);
t144 = cos(qJ(3));
t153 = t135 * t144;
t136 = cos(t137);
t152 = t136 * t144;
t139 = sin(pkin(6));
t142 = sin(qJ(3));
t151 = t139 * t142;
t150 = t139 * t144;
t145 = cos(qJ(2));
t149 = t139 * t145;
t141 = cos(pkin(6));
t143 = sin(qJ(2));
t148 = t141 * t143;
t147 = t141 * t145;
t146 = t144 * t145;
t140 = cos(pkin(11));
t138 = sin(pkin(11));
t133 = t141 * t142 + t143 * t150;
t132 = t141 * t144 - t143 * t151;
t131 = -t138 * t148 + t140 * t145;
t130 = t138 * t147 + t140 * t143;
t129 = t138 * t145 + t140 * t148;
t128 = t138 * t143 - t140 * t147;
t127 = t131 * t144 + t138 * t151;
t126 = -t131 * t142 + t138 * t150;
t125 = t129 * t144 - t140 * t151;
t124 = -t129 * t142 - t140 * t150;
t123 = -t133 * t136 + t135 * t149;
t122 = -t133 * t135 - t136 * t149;
t121 = -t127 * t136 - t130 * t135;
t120 = -t127 * t135 + t130 * t136;
t119 = -t125 * t136 - t128 * t135;
t118 = -t125 * t135 + t128 * t136;
t1 = [0, -t130 * t152 + t131 * t135, t126 * t136, 0, t120, t120; 0, -t128 * t152 + t129 * t135, t124 * t136, 0, t118, t118; 0 (t135 * t143 + t136 * t146) * t139, t132 * t136, 0, t122, t122; 0, t130 * t153 + t131 * t136, -t126 * t135, 0, t121, t121; 0, t128 * t153 + t129 * t136, -t124 * t135, 0, t119, t119; 0 (-t135 * t146 + t136 * t143) * t139, -t132 * t135, 0, t123, t123; 0, -t130 * t142, t127, 0, 0, 0; 0, -t128 * t142, t125, 0, 0, 0; 0, t142 * t149, t133, 0, 0, 0;];
JR_rot  = t1;
