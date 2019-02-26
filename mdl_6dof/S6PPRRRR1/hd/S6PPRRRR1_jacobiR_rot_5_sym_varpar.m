% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:43
% EndTime: 2019-02-26 19:42:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (135->27), mult. (307->58), div. (0->0), fcn. (430->12), ass. (0->41)
t144 = sin(pkin(12));
t150 = cos(pkin(6));
t160 = t144 * t150;
t145 = sin(pkin(7));
t146 = sin(pkin(6));
t159 = t145 * t146;
t158 = t145 * t150;
t149 = cos(pkin(7));
t157 = t146 * t149;
t147 = cos(pkin(13));
t156 = t147 * t149;
t148 = cos(pkin(12));
t155 = t148 * t150;
t143 = sin(pkin(13));
t136 = -t144 * t143 + t147 * t155;
t154 = t136 * t149 - t148 * t159;
t138 = -t148 * t143 - t147 * t160;
t153 = t138 * t149 + t144 * t159;
t152 = cos(qJ(3));
t151 = sin(qJ(3));
t142 = qJ(4) + qJ(5);
t141 = cos(t142);
t140 = sin(t142);
t139 = -t143 * t160 + t148 * t147;
t137 = t143 * t155 + t144 * t147;
t135 = -t147 * t159 + t150 * t149;
t134 = -t138 * t145 + t144 * t157;
t133 = -t136 * t145 - t148 * t157;
t132 = t151 * t158 + (t143 * t152 + t151 * t156) * t146;
t131 = t152 * t158 + (-t143 * t151 + t152 * t156) * t146;
t130 = t139 * t152 + t153 * t151;
t129 = -t139 * t151 + t153 * t152;
t128 = t137 * t152 + t154 * t151;
t127 = -t137 * t151 + t154 * t152;
t126 = -t132 * t141 - t135 * t140;
t125 = -t132 * t140 + t135 * t141;
t124 = -t130 * t141 - t134 * t140;
t123 = -t130 * t140 + t134 * t141;
t122 = -t128 * t141 - t133 * t140;
t121 = -t128 * t140 + t133 * t141;
t1 = [0, 0, t129 * t141, t123, t123, 0; 0, 0, t127 * t141, t121, t121, 0; 0, 0, t131 * t141, t125, t125, 0; 0, 0, -t129 * t140, t124, t124, 0; 0, 0, -t127 * t140, t122, t122, 0; 0, 0, -t131 * t140, t126, t126, 0; 0, 0, t130, 0, 0, 0; 0, 0, t128, 0, 0, 0; 0, 0, t132, 0, 0, 0;];
JR_rot  = t1;
