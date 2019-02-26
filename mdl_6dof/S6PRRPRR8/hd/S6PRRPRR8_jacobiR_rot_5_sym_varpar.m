% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR8_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:07
% EndTime: 2019-02-26 20:08:07
% DurationCPUTime: 0.12s
% Computational Cost: add. (97->38), mult. (304->88), div. (0->0), fcn. (425->12), ass. (0->45)
t153 = sin(pkin(7));
t154 = sin(pkin(6));
t179 = t153 * t154;
t157 = cos(pkin(6));
t178 = t153 * t157;
t158 = sin(qJ(5));
t177 = t153 * t158;
t161 = cos(qJ(5));
t176 = t153 * t161;
t156 = cos(pkin(7));
t175 = t154 * t156;
t159 = sin(qJ(3));
t174 = t156 * t159;
t162 = cos(qJ(3));
t173 = t156 * t162;
t160 = sin(qJ(2));
t172 = t157 * t160;
t163 = cos(qJ(2));
t171 = t157 * t163;
t170 = t159 * t160;
t169 = t159 * t163;
t168 = t160 * t162;
t167 = t162 * t163;
t166 = t160 * t179;
t152 = sin(pkin(12));
t155 = cos(pkin(12));
t147 = -t152 * t160 + t155 * t171;
t165 = -t147 * t156 + t155 * t179;
t149 = -t152 * t171 - t155 * t160;
t164 = t149 * t156 + t152 * t179;
t150 = -t152 * t172 + t155 * t163;
t148 = t152 * t163 + t155 * t172;
t146 = t157 * t156 - t163 * t179;
t145 = (t156 * t168 + t169) * t154;
t144 = -t149 * t153 + t152 * t175;
t143 = -t147 * t153 - t155 * t175;
t142 = t159 * t178 + (t156 * t169 + t168) * t154;
t141 = -t162 * t178 + (-t156 * t167 + t170) * t154;
t140 = t149 * t159 + t150 * t173;
t139 = t147 * t159 + t148 * t173;
t138 = t150 * t162 + t164 * t159;
t137 = t150 * t159 - t164 * t162;
t136 = t148 * t162 - t165 * t159;
t135 = t148 * t159 + t165 * t162;
t1 = [0, t140 * t158 + t150 * t176, t138 * t158, 0, t137 * t161 - t144 * t158, 0; 0, t139 * t158 + t148 * t176, t136 * t158, 0, t135 * t161 - t143 * t158, 0; 0, t145 * t158 + t161 * t166, t142 * t158, 0, t141 * t161 - t146 * t158, 0; 0, t140 * t161 - t150 * t177, t138 * t161, 0, -t137 * t158 - t144 * t161, 0; 0, t139 * t161 - t148 * t177, t136 * t161, 0, -t135 * t158 - t143 * t161, 0; 0, t145 * t161 - t158 * t166, t142 * t161, 0, -t141 * t158 - t146 * t161, 0; 0, t149 * t162 - t150 * t174, -t137, 0, 0, 0; 0, t147 * t162 - t148 * t174, -t135, 0, 0, 0; 0 (-t156 * t170 + t167) * t154, -t141, 0, 0, 0;];
JR_rot  = t1;
