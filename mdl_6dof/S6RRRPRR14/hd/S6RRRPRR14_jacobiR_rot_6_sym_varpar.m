% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR14_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:40
% EndTime: 2019-02-26 22:23:40
% DurationCPUTime: 0.10s
% Computational Cost: add. (141->33), mult. (275->62), div. (0->0), fcn. (402->10), ass. (0->39)
t153 = cos(pkin(6));
t155 = sin(qJ(2));
t159 = cos(qJ(1));
t162 = t159 * t155;
t156 = sin(qJ(1));
t158 = cos(qJ(2));
t163 = t156 * t158;
t145 = t153 * t162 + t163;
t154 = sin(qJ(3));
t157 = cos(qJ(3));
t152 = sin(pkin(6));
t166 = t152 * t159;
t137 = t145 * t154 + t157 * t166;
t161 = t159 * t158;
t164 = t156 * t155;
t144 = -t153 * t161 + t164;
t151 = qJ(5) + qJ(6);
t149 = sin(t151);
t150 = cos(t151);
t132 = -t137 * t149 - t144 * t150;
t131 = t137 * t150 - t144 * t149;
t171 = t149 * t154;
t170 = t150 * t154;
t169 = t152 * t154;
t168 = t152 * t157;
t167 = t152 * t158;
t165 = t154 * t158;
t160 = -t145 * t157 + t154 * t166;
t147 = -t153 * t164 + t161;
t146 = t153 * t163 + t162;
t143 = t153 * t154 + t155 * t168;
t142 = -t153 * t157 + t155 * t169;
t141 = t147 * t157 + t156 * t169;
t140 = t147 * t154 - t156 * t168;
t136 = -t142 * t149 + t150 * t167;
t135 = t142 * t150 + t149 * t167;
t134 = t140 * t149 + t146 * t150;
t133 = t140 * t150 - t146 * t149;
t1 = [t132, -t146 * t171 + t147 * t150, t141 * t149, 0, t133, t133; t134, -t144 * t171 + t145 * t150, -t160 * t149, 0, t131, t131; 0 (t149 * t165 + t150 * t155) * t152, t143 * t149, 0, t135, t135; -t131, -t146 * t170 - t147 * t149, t141 * t150, 0, -t134, -t134; t133, -t144 * t170 - t145 * t149, -t160 * t150, 0, t132, t132; 0 (-t149 * t155 + t150 * t165) * t152, t143 * t150, 0, t136, t136; t160, -t146 * t157, -t140, 0, 0, 0; t141, -t144 * t157, -t137, 0, 0, 0; 0, t157 * t167, -t142, 0, 0, 0;];
JR_rot  = t1;
