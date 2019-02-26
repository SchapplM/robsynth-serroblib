% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR3
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
% Datum: 2019-02-26 19:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:54:49
% EndTime: 2019-02-26 19:54:49
% DurationCPUTime: 0.08s
% Computational Cost: add. (186->31), mult. (214->65), div. (0->0), fcn. (313->10), ass. (0->37)
t146 = sin(pkin(11));
t148 = cos(pkin(11));
t153 = cos(qJ(2));
t149 = cos(pkin(6));
t151 = sin(qJ(2));
t157 = t149 * t151;
t139 = t146 * t153 + t148 * t157;
t145 = pkin(12) + qJ(4) + qJ(5);
t143 = sin(t145);
t144 = cos(t145);
t147 = sin(pkin(6));
t159 = t147 * t148;
t131 = -t139 * t143 - t144 * t159;
t150 = sin(qJ(6));
t165 = t131 * t150;
t141 = -t146 * t157 + t148 * t153;
t160 = t146 * t147;
t133 = -t141 * t143 + t144 * t160;
t164 = t133 * t150;
t158 = t147 * t151;
t136 = -t143 * t158 + t149 * t144;
t163 = t136 * t150;
t162 = t144 * t150;
t152 = cos(qJ(6));
t161 = t144 * t152;
t156 = t149 * t153;
t155 = t150 * t153;
t154 = t152 * t153;
t140 = t146 * t156 + t148 * t151;
t138 = t146 * t151 - t148 * t156;
t137 = t149 * t143 + t144 * t158;
t135 = t136 * t152;
t134 = t141 * t144 + t143 * t160;
t132 = t139 * t144 - t143 * t159;
t130 = t133 * t152;
t129 = t131 * t152;
t1 = [0, -t140 * t161 + t141 * t150, 0, t130, t130, -t134 * t150 + t140 * t152; 0, -t138 * t161 + t139 * t150, 0, t129, t129, -t132 * t150 + t138 * t152; 0 (t144 * t154 + t150 * t151) * t147, 0, t135, t135, -t137 * t150 - t147 * t154; 0, t140 * t162 + t141 * t152, 0, -t164, -t164, -t134 * t152 - t140 * t150; 0, t138 * t162 + t139 * t152, 0, -t165, -t165, -t132 * t152 - t138 * t150; 0 (-t144 * t155 + t151 * t152) * t147, 0, -t163, -t163, -t137 * t152 + t147 * t155; 0, -t140 * t143, 0, t134, t134, 0; 0, -t138 * t143, 0, t132, t132, 0; 0, t147 * t153 * t143, 0, t137, t137, 0;];
JR_rot  = t1;
