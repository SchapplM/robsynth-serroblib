% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:18
% EndTime: 2019-02-26 19:50:18
% DurationCPUTime: 0.13s
% Computational Cost: add. (114->28), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->33)
t155 = sin(pkin(11));
t158 = cos(pkin(11));
t163 = sin(qJ(2));
t166 = cos(qJ(2));
t151 = t163 * t155 - t166 * t158;
t157 = sin(pkin(6));
t162 = sin(qJ(4));
t176 = t157 * t162;
t165 = cos(qJ(4));
t175 = t157 * t165;
t161 = sin(qJ(5));
t174 = t161 * t165;
t164 = cos(qJ(5));
t172 = t164 * t165;
t160 = cos(pkin(6));
t168 = t166 * t155 + t163 * t158;
t150 = t168 * t160;
t156 = sin(pkin(10));
t159 = cos(pkin(10));
t170 = t159 * t150 - t156 * t151;
t169 = -t156 * t150 - t159 * t151;
t167 = t151 * t160;
t149 = t168 * t157;
t148 = t151 * t157;
t146 = t149 * t165 + t160 * t162;
t145 = -t149 * t162 + t160 * t165;
t143 = t156 * t167 - t159 * t168;
t140 = -t156 * t168 - t159 * t167;
t138 = t156 * t176 + t165 * t169;
t137 = t156 * t175 - t162 * t169;
t136 = -t159 * t176 + t165 * t170;
t135 = -t159 * t175 - t162 * t170;
t1 = [0, t143 * t172 + t161 * t169, 0, t137 * t164, -t138 * t161 - t143 * t164, 0; 0, t140 * t172 + t161 * t170, 0, t135 * t164, -t136 * t161 - t140 * t164, 0; 0, -t148 * t172 + t149 * t161, 0, t145 * t164, -t146 * t161 + t148 * t164, 0; 0, -t143 * t174 + t164 * t169, 0, -t137 * t161, -t138 * t164 + t143 * t161, 0; 0, -t140 * t174 + t164 * t170, 0, -t135 * t161, -t136 * t164 + t140 * t161, 0; 0, t148 * t174 + t149 * t164, 0, -t145 * t161, -t146 * t164 - t148 * t161, 0; 0, t143 * t162, 0, t138, 0, 0; 0, t140 * t162, 0, t136, 0, 0; 0, -t148 * t162, 0, t146, 0, 0;];
JR_rot  = t1;
