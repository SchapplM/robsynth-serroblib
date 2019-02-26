% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRP3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:02:23
% EndTime: 2019-02-26 20:02:23
% DurationCPUTime: 0.10s
% Computational Cost: add. (81->25), mult. (163->63), div. (0->0), fcn. (238->10), ass. (0->30)
t144 = pkin(11) + qJ(5);
t142 = sin(t144);
t151 = cos(qJ(3));
t160 = t142 * t151;
t143 = cos(t144);
t159 = t143 * t151;
t146 = sin(pkin(6));
t149 = sin(qJ(3));
t158 = t146 * t149;
t157 = t146 * t151;
t152 = cos(qJ(2));
t156 = t146 * t152;
t148 = cos(pkin(6));
t150 = sin(qJ(2));
t155 = t148 * t150;
t154 = t148 * t152;
t153 = t151 * t152;
t147 = cos(pkin(10));
t145 = sin(pkin(10));
t140 = t148 * t149 + t150 * t157;
t139 = t148 * t151 - t150 * t158;
t138 = -t145 * t155 + t147 * t152;
t137 = t145 * t154 + t147 * t150;
t136 = t145 * t152 + t147 * t155;
t135 = t145 * t150 - t147 * t154;
t134 = t138 * t151 + t145 * t158;
t133 = -t138 * t149 + t145 * t157;
t132 = t136 * t151 - t147 * t158;
t131 = -t136 * t149 - t147 * t157;
t1 = [0, -t137 * t159 + t138 * t142, t133 * t143, 0, -t134 * t142 + t137 * t143, 0; 0, -t135 * t159 + t136 * t142, t131 * t143, 0, -t132 * t142 + t135 * t143, 0; 0 (t142 * t150 + t143 * t153) * t146, t139 * t143, 0, -t140 * t142 - t143 * t156, 0; 0, -t137 * t149, t134, 0, 0, 0; 0, -t135 * t149, t132, 0, 0, 0; 0, t149 * t156, t140, 0, 0, 0; 0, -t137 * t160 - t138 * t143, t133 * t142, 0, t134 * t143 + t137 * t142, 0; 0, -t135 * t160 - t136 * t143, t131 * t142, 0, t132 * t143 + t135 * t142, 0; 0 (t142 * t153 - t143 * t150) * t146, t139 * t142, 0, t140 * t143 - t142 * t156, 0;];
JR_rot  = t1;
