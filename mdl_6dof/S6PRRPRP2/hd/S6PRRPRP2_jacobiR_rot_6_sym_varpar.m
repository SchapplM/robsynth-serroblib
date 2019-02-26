% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP2
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
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRP2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:48
% EndTime: 2019-02-26 20:01:48
% DurationCPUTime: 0.07s
% Computational Cost: add. (90->25), mult. (163->65), div. (0->0), fcn. (238->10), ass. (0->31)
t146 = qJ(3) + pkin(11);
t145 = cos(t146);
t151 = sin(qJ(5));
t163 = t145 * t151;
t153 = cos(qJ(5));
t162 = t145 * t153;
t147 = sin(pkin(10));
t148 = sin(pkin(6));
t161 = t147 * t148;
t149 = cos(pkin(10));
t160 = t148 * t149;
t152 = sin(qJ(2));
t159 = t148 * t152;
t150 = cos(pkin(6));
t158 = t150 * t152;
t154 = cos(qJ(2));
t157 = t150 * t154;
t156 = t151 * t154;
t155 = t153 * t154;
t144 = sin(t146);
t142 = -t147 * t158 + t149 * t154;
t141 = t147 * t157 + t149 * t152;
t140 = t147 * t154 + t149 * t158;
t139 = t147 * t152 - t149 * t157;
t138 = t150 * t144 + t145 * t159;
t137 = -t144 * t159 + t150 * t145;
t136 = t142 * t145 + t144 * t161;
t135 = -t142 * t144 + t145 * t161;
t134 = t140 * t145 - t144 * t160;
t133 = -t140 * t144 - t145 * t160;
t1 = [0, -t141 * t162 + t142 * t151, t135 * t153, 0, -t136 * t151 + t141 * t153, 0; 0, -t139 * t162 + t140 * t151, t133 * t153, 0, -t134 * t151 + t139 * t153, 0; 0 (t145 * t155 + t151 * t152) * t148, t137 * t153, 0, -t138 * t151 - t148 * t155, 0; 0, -t141 * t144, t136, 0, 0, 0; 0, -t139 * t144, t134, 0, 0, 0; 0, t148 * t154 * t144, t138, 0, 0, 0; 0, -t141 * t163 - t142 * t153, t135 * t151, 0, t136 * t153 + t141 * t151, 0; 0, -t139 * t163 - t140 * t153, t133 * t151, 0, t134 * t153 + t139 * t151, 0; 0 (t145 * t156 - t152 * t153) * t148, t137 * t151, 0, t138 * t153 - t148 * t156, 0;];
JR_rot  = t1;
