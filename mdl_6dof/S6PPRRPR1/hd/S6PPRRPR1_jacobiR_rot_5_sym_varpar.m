% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRPR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:29
% EndTime: 2019-02-26 19:40:29
% DurationCPUTime: 0.10s
% Computational Cost: add. (113->32), mult. (338->75), div. (0->0), fcn. (466->14), ass. (0->41)
t152 = sin(pkin(13));
t164 = cos(qJ(4));
t175 = t152 * t164;
t154 = sin(pkin(11));
t161 = cos(pkin(6));
t174 = t154 * t161;
t155 = sin(pkin(7));
t156 = sin(pkin(6));
t173 = t155 * t156;
t172 = t155 * t161;
t160 = cos(pkin(7));
t171 = t156 * t160;
t157 = cos(pkin(13));
t170 = t157 * t164;
t158 = cos(pkin(12));
t169 = t158 * t160;
t159 = cos(pkin(11));
t168 = t159 * t161;
t153 = sin(pkin(12));
t148 = -t153 * t154 + t158 * t168;
t167 = t148 * t160 - t159 * t173;
t150 = -t153 * t159 - t158 * t174;
t166 = t150 * t160 + t154 * t173;
t165 = cos(qJ(3));
t163 = sin(qJ(3));
t162 = sin(qJ(4));
t151 = -t153 * t174 + t158 * t159;
t149 = t153 * t168 + t154 * t158;
t147 = -t158 * t173 + t160 * t161;
t146 = -t150 * t155 + t154 * t171;
t145 = -t148 * t155 - t159 * t171;
t144 = t163 * t172 + (t153 * t165 + t163 * t169) * t156;
t143 = t165 * t172 + (-t153 * t163 + t165 * t169) * t156;
t142 = -t144 * t162 + t147 * t164;
t141 = t151 * t165 + t166 * t163;
t140 = -t151 * t163 + t166 * t165;
t139 = t149 * t165 + t167 * t163;
t138 = -t149 * t163 + t167 * t165;
t137 = -t141 * t162 + t146 * t164;
t136 = -t139 * t162 + t145 * t164;
t1 = [0, 0, t140 * t170 + t141 * t152, t137 * t157, 0, 0; 0, 0, t138 * t170 + t139 * t152, t136 * t157, 0, 0; 0, 0, t143 * t170 + t144 * t152, t142 * t157, 0, 0; 0, 0, -t140 * t175 + t141 * t157, -t137 * t152, 0, 0; 0, 0, -t138 * t175 + t139 * t157, -t136 * t152, 0, 0; 0, 0, -t143 * t175 + t144 * t157, -t142 * t152, 0, 0; 0, 0, t140 * t162, t141 * t164 + t146 * t162, 0, 0; 0, 0, t138 * t162, t139 * t164 + t145 * t162, 0, 0; 0, 0, t143 * t162, t144 * t164 + t147 * t162, 0, 0;];
JR_rot  = t1;
