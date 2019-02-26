% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP1
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
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:13
% EndTime: 2019-02-26 20:15:13
% DurationCPUTime: 0.13s
% Computational Cost: add. (129->31), mult. (214->65), div. (0->0), fcn. (313->10), ass. (0->37)
t154 = sin(pkin(11));
t156 = cos(pkin(11));
t161 = cos(qJ(2));
t157 = cos(pkin(6));
t159 = sin(qJ(2));
t165 = t157 * t159;
t147 = t154 * t161 + t156 * t165;
t153 = qJ(3) + qJ(4);
t151 = sin(t153);
t152 = cos(t153);
t155 = sin(pkin(6));
t167 = t155 * t156;
t139 = -t147 * t151 - t152 * t167;
t158 = sin(qJ(5));
t173 = t139 * t158;
t149 = -t154 * t165 + t156 * t161;
t168 = t154 * t155;
t141 = -t149 * t151 + t152 * t168;
t172 = t141 * t158;
t166 = t155 * t159;
t144 = -t151 * t166 + t152 * t157;
t171 = t144 * t158;
t170 = t152 * t158;
t160 = cos(qJ(5));
t169 = t152 * t160;
t164 = t157 * t161;
t163 = t158 * t161;
t162 = t160 * t161;
t148 = t154 * t164 + t156 * t159;
t146 = t154 * t159 - t156 * t164;
t145 = t151 * t157 + t152 * t166;
t143 = t144 * t160;
t142 = t149 * t152 + t151 * t168;
t140 = t147 * t152 - t151 * t167;
t138 = t141 * t160;
t137 = t139 * t160;
t1 = [0, -t148 * t169 + t149 * t158, t138, t138, -t142 * t158 + t148 * t160, 0; 0, -t146 * t169 + t147 * t158, t137, t137, -t140 * t158 + t146 * t160, 0; 0 (t152 * t162 + t158 * t159) * t155, t143, t143, -t145 * t158 - t155 * t162, 0; 0, t148 * t170 + t149 * t160, -t172, -t172, -t142 * t160 - t148 * t158, 0; 0, t146 * t170 + t147 * t160, -t173, -t173, -t140 * t160 - t146 * t158, 0; 0 (-t152 * t163 + t159 * t160) * t155, -t171, -t171, -t145 * t160 + t155 * t163, 0; 0, -t148 * t151, t142, t142, 0, 0; 0, -t146 * t151, t140, t140, 0, 0; 0, t155 * t161 * t151, t145, t145, 0, 0;];
JR_rot  = t1;
