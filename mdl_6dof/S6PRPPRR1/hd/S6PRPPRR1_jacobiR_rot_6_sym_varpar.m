% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPPRR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (153->29), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->34)
t153 = sin(pkin(11));
t156 = cos(pkin(11));
t160 = sin(qJ(2));
t162 = cos(qJ(2));
t146 = t160 * t153 - t162 * t156;
t152 = pkin(12) + qJ(5);
t151 = cos(t152);
t159 = sin(qJ(6));
t172 = t151 * t159;
t161 = cos(qJ(6));
t171 = t151 * t161;
t154 = sin(pkin(10));
t155 = sin(pkin(6));
t170 = t154 * t155;
t157 = cos(pkin(10));
t169 = t155 * t157;
t158 = cos(pkin(6));
t164 = t162 * t153 + t160 * t156;
t145 = t164 * t158;
t166 = t157 * t145 - t154 * t146;
t165 = -t154 * t145 - t157 * t146;
t163 = t146 * t158;
t150 = sin(t152);
t144 = t164 * t155;
t143 = t146 * t155;
t141 = t144 * t151 + t158 * t150;
t140 = -t144 * t150 + t158 * t151;
t138 = t154 * t163 - t157 * t164;
t135 = -t154 * t164 - t157 * t163;
t133 = t150 * t170 + t151 * t165;
t132 = -t150 * t165 + t151 * t170;
t131 = -t150 * t169 + t151 * t166;
t130 = -t150 * t166 - t151 * t169;
t1 = [0, t138 * t171 + t159 * t165, 0, 0, t132 * t161, -t133 * t159 - t138 * t161; 0, t135 * t171 + t159 * t166, 0, 0, t130 * t161, -t131 * t159 - t135 * t161; 0, -t143 * t171 + t144 * t159, 0, 0, t140 * t161, -t141 * t159 + t143 * t161; 0, -t138 * t172 + t161 * t165, 0, 0, -t132 * t159, -t133 * t161 + t138 * t159; 0, -t135 * t172 + t161 * t166, 0, 0, -t130 * t159, -t131 * t161 + t135 * t159; 0, t143 * t172 + t144 * t161, 0, 0, -t140 * t159, -t141 * t161 - t143 * t159; 0, t138 * t150, 0, 0, t133, 0; 0, t135 * t150, 0, 0, t131, 0; 0, -t143 * t150, 0, 0, t141, 0;];
JR_rot  = t1;
