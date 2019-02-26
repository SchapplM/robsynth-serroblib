% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:07
% EndTime: 2019-02-26 20:11:07
% DurationCPUTime: 0.12s
% Computational Cost: add. (165->32), mult. (214->64), div. (0->0), fcn. (313->10), ass. (0->38)
t155 = sin(pkin(11));
t157 = cos(pkin(11));
t160 = cos(qJ(2));
t158 = cos(pkin(6));
t159 = sin(qJ(2));
t162 = t158 * t159;
t145 = t155 * t160 + t157 * t162;
t154 = qJ(3) + qJ(4);
t151 = sin(t154);
t152 = cos(t154);
t156 = sin(pkin(6));
t165 = t156 * t157;
t137 = -t145 * t151 - t152 * t165;
t153 = pkin(12) + qJ(6);
t149 = sin(t153);
t172 = t137 * t149;
t147 = -t155 * t162 + t157 * t160;
t166 = t155 * t156;
t139 = -t147 * t151 + t152 * t166;
t171 = t139 * t149;
t164 = t156 * t159;
t142 = -t151 * t164 + t158 * t152;
t170 = t142 * t149;
t169 = t149 * t152;
t150 = cos(t153);
t168 = t150 * t152;
t167 = t152 * t160;
t163 = t156 * t160;
t161 = t158 * t160;
t146 = t155 * t161 + t157 * t159;
t144 = t155 * t159 - t157 * t161;
t143 = t158 * t151 + t152 * t164;
t141 = t142 * t150;
t140 = t147 * t152 + t151 * t166;
t138 = t145 * t152 - t151 * t165;
t136 = t139 * t150;
t135 = t137 * t150;
t1 = [0, -t146 * t168 + t147 * t149, t136, t136, 0, -t140 * t149 + t146 * t150; 0, -t144 * t168 + t145 * t149, t135, t135, 0, -t138 * t149 + t144 * t150; 0 (t149 * t159 + t150 * t167) * t156, t141, t141, 0, -t143 * t149 - t150 * t163; 0, t146 * t169 + t147 * t150, -t171, -t171, 0, -t140 * t150 - t146 * t149; 0, t144 * t169 + t145 * t150, -t172, -t172, 0, -t138 * t150 - t144 * t149; 0 (-t149 * t167 + t150 * t159) * t156, -t170, -t170, 0, -t143 * t150 + t149 * t163; 0, -t146 * t151, t140, t140, 0, 0; 0, -t144 * t151, t138, t138, 0, 0; 0, t151 * t163, t143, t143, 0, 0;];
JR_rot  = t1;
