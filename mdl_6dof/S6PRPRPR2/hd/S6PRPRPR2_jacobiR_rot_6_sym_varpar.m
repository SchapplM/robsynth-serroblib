% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:03
% EndTime: 2019-02-26 19:47:03
% DurationCPUTime: 0.14s
% Computational Cost: add. (144->29), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->34)
t155 = sin(pkin(11));
t158 = cos(pkin(11));
t162 = sin(qJ(2));
t164 = cos(qJ(2));
t148 = t162 * t155 - t164 * t158;
t154 = pkin(12) + qJ(6);
t152 = sin(t154);
t163 = cos(qJ(4));
t174 = t152 * t163;
t153 = cos(t154);
t173 = t153 * t163;
t157 = sin(pkin(6));
t161 = sin(qJ(4));
t172 = t157 * t161;
t171 = t157 * t163;
t160 = cos(pkin(6));
t166 = t164 * t155 + t162 * t158;
t147 = t166 * t160;
t156 = sin(pkin(10));
t159 = cos(pkin(10));
t168 = t159 * t147 - t156 * t148;
t167 = -t156 * t147 - t159 * t148;
t165 = t148 * t160;
t146 = t166 * t157;
t145 = t148 * t157;
t143 = t146 * t163 + t160 * t161;
t142 = -t146 * t161 + t160 * t163;
t140 = t156 * t165 - t159 * t166;
t137 = -t156 * t166 - t159 * t165;
t135 = t156 * t172 + t163 * t167;
t134 = t156 * t171 - t161 * t167;
t133 = -t159 * t172 + t163 * t168;
t132 = -t159 * t171 - t161 * t168;
t1 = [0, t140 * t173 + t152 * t167, 0, t134 * t153, 0, -t135 * t152 - t140 * t153; 0, t137 * t173 + t152 * t168, 0, t132 * t153, 0, -t133 * t152 - t137 * t153; 0, -t145 * t173 + t146 * t152, 0, t142 * t153, 0, -t143 * t152 + t145 * t153; 0, -t140 * t174 + t153 * t167, 0, -t134 * t152, 0, -t135 * t153 + t140 * t152; 0, -t137 * t174 + t153 * t168, 0, -t132 * t152, 0, -t133 * t153 + t137 * t152; 0, t145 * t174 + t146 * t153, 0, -t142 * t152, 0, -t143 * t153 - t145 * t152; 0, t140 * t161, 0, t135, 0, 0; 0, t137 * t161, 0, t133, 0, 0; 0, -t145 * t161, 0, t143, 0, 0;];
JR_rot  = t1;
