% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRP3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:21
% EndTime: 2019-02-26 20:16:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (116->28), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->36)
t140 = qJ(4) + qJ(5);
t138 = sin(t140);
t147 = cos(qJ(3));
t156 = t138 * t147;
t139 = cos(t140);
t155 = t139 * t147;
t142 = sin(pkin(6));
t145 = sin(qJ(3));
t154 = t142 * t145;
t153 = t142 * t147;
t148 = cos(qJ(2));
t152 = t142 * t148;
t144 = cos(pkin(6));
t146 = sin(qJ(2));
t151 = t144 * t146;
t150 = t144 * t148;
t149 = t147 * t148;
t143 = cos(pkin(11));
t141 = sin(pkin(11));
t136 = t144 * t145 + t146 * t153;
t135 = t144 * t147 - t146 * t154;
t134 = -t141 * t151 + t143 * t148;
t133 = t141 * t150 + t143 * t146;
t132 = t141 * t148 + t143 * t151;
t131 = t141 * t146 - t143 * t150;
t130 = t134 * t147 + t141 * t154;
t129 = -t134 * t145 + t141 * t153;
t128 = t132 * t147 - t143 * t154;
t127 = -t132 * t145 - t143 * t153;
t126 = -t136 * t139 + t138 * t152;
t125 = -t136 * t138 - t139 * t152;
t124 = -t130 * t139 - t133 * t138;
t123 = -t130 * t138 + t133 * t139;
t122 = -t128 * t139 - t131 * t138;
t121 = -t128 * t138 + t131 * t139;
t1 = [0, -t133 * t155 + t134 * t138, t129 * t139, t123, t123, 0; 0, -t131 * t155 + t132 * t138, t127 * t139, t121, t121, 0; 0 (t138 * t146 + t139 * t149) * t142, t135 * t139, t125, t125, 0; 0, t133 * t156 + t134 * t139, -t129 * t138, t124, t124, 0; 0, t131 * t156 + t132 * t139, -t127 * t138, t122, t122, 0; 0 (-t138 * t149 + t139 * t146) * t142, -t135 * t138, t126, t126, 0; 0, -t133 * t145, t130, 0, 0, 0; 0, -t131 * t145, t128, 0, 0, 0; 0, t145 * t152, t136, 0, 0, 0;];
JR_rot  = t1;
