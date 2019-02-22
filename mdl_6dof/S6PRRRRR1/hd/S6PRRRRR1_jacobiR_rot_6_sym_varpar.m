% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR1
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
% Datum: 2019-02-22 10:02
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:02:00
% EndTime: 2019-02-22 10:02:01
% DurationCPUTime: 0.15s
% Computational Cost: add. (240->34), mult. (265->65), div. (0->0), fcn. (388->10), ass. (0->37)
t152 = sin(pkin(12));
t154 = cos(pkin(12));
t159 = cos(qJ(2));
t155 = cos(pkin(6));
t157 = sin(qJ(2));
t163 = t155 * t157;
t145 = t152 * t159 + t154 * t163;
t151 = qJ(3) + qJ(4) + qJ(5);
t149 = sin(t151);
t150 = cos(t151);
t153 = sin(pkin(6));
t165 = t153 * t154;
t137 = -t145 * t149 - t150 * t165;
t156 = sin(qJ(6));
t171 = t137 * t156;
t147 = -t152 * t163 + t154 * t159;
t166 = t152 * t153;
t139 = -t147 * t149 + t150 * t166;
t170 = t139 * t156;
t164 = t153 * t157;
t142 = -t149 * t164 + t155 * t150;
t169 = t142 * t156;
t168 = t150 * t156;
t158 = cos(qJ(6));
t167 = t150 * t158;
t162 = t155 * t159;
t161 = t156 * t159;
t160 = t158 * t159;
t146 = t152 * t162 + t154 * t157;
t144 = t152 * t157 - t154 * t162;
t143 = t155 * t149 + t150 * t164;
t141 = t142 * t158;
t140 = t147 * t150 + t149 * t166;
t138 = t145 * t150 - t149 * t165;
t136 = t139 * t158;
t135 = t137 * t158;
t1 = [0, -t146 * t167 + t147 * t156, t136, t136, t136, -t140 * t156 + t146 * t158; 0, -t144 * t167 + t145 * t156, t135, t135, t135, -t138 * t156 + t144 * t158; 0 (t150 * t160 + t156 * t157) * t153, t141, t141, t141, -t143 * t156 - t153 * t160; 0, t146 * t168 + t147 * t158, -t170, -t170, -t170, -t140 * t158 - t146 * t156; 0, t144 * t168 + t145 * t158, -t171, -t171, -t171, -t138 * t158 - t144 * t156; 0 (-t150 * t161 + t157 * t158) * t153, -t169, -t169, -t169, -t143 * t158 + t153 * t161; 0, -t146 * t149, t140, t140, t140, 0; 0, -t144 * t149, t138, t138, t138, 0; 0, t153 * t159 * t149, t143, t143, t143, 0;];
JR_rot  = t1;
