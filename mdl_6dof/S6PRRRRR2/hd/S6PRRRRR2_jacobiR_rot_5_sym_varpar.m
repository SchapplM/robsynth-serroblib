% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR2
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

function JR_rot = S6PRRRRR2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:15
% EndTime: 2019-02-26 20:19:15
% DurationCPUTime: 0.08s
% Computational Cost: add. (129->31), mult. (214->65), div. (0->0), fcn. (313->10), ass. (0->37)
t145 = sin(pkin(12));
t147 = cos(pkin(12));
t152 = cos(qJ(2));
t148 = cos(pkin(6));
t150 = sin(qJ(2));
t156 = t148 * t150;
t138 = t145 * t152 + t147 * t156;
t144 = qJ(3) + qJ(4);
t142 = sin(t144);
t143 = cos(t144);
t146 = sin(pkin(6));
t158 = t146 * t147;
t130 = -t138 * t142 - t143 * t158;
t149 = sin(qJ(5));
t164 = t130 * t149;
t140 = -t145 * t156 + t147 * t152;
t159 = t145 * t146;
t132 = -t140 * t142 + t143 * t159;
t163 = t132 * t149;
t157 = t146 * t150;
t135 = -t142 * t157 + t148 * t143;
t162 = t135 * t149;
t161 = t143 * t149;
t151 = cos(qJ(5));
t160 = t143 * t151;
t155 = t148 * t152;
t154 = t149 * t152;
t153 = t151 * t152;
t139 = t145 * t155 + t147 * t150;
t137 = t145 * t150 - t147 * t155;
t136 = t148 * t142 + t143 * t157;
t134 = t135 * t151;
t133 = t140 * t143 + t142 * t159;
t131 = t138 * t143 - t142 * t158;
t129 = t132 * t151;
t128 = t130 * t151;
t1 = [0, -t139 * t160 + t140 * t149, t129, t129, -t133 * t149 + t139 * t151, 0; 0, -t137 * t160 + t138 * t149, t128, t128, -t131 * t149 + t137 * t151, 0; 0 (t143 * t153 + t149 * t150) * t146, t134, t134, -t136 * t149 - t146 * t153, 0; 0, t139 * t161 + t140 * t151, -t163, -t163, -t133 * t151 - t139 * t149, 0; 0, t137 * t161 + t138 * t151, -t164, -t164, -t131 * t151 - t137 * t149, 0; 0 (-t143 * t154 + t150 * t151) * t146, -t162, -t162, -t136 * t151 + t146 * t154, 0; 0, -t139 * t142, t133, t133, 0, 0; 0, -t137 * t142, t131, t131, 0, 0; 0, t146 * t152 * t142, t136, t136, 0, 0;];
JR_rot  = t1;
