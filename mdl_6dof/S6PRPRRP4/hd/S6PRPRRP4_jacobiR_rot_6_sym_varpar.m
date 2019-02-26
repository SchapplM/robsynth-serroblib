% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP4
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
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRP4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:01
% EndTime: 2019-02-26 19:52:01
% DurationCPUTime: 0.12s
% Computational Cost: add. (90->25), mult. (163->65), div. (0->0), fcn. (238->10), ass. (0->31)
t143 = pkin(11) + qJ(4);
t142 = cos(t143);
t148 = sin(qJ(5));
t160 = t142 * t148;
t150 = cos(qJ(5));
t159 = t142 * t150;
t144 = sin(pkin(10));
t145 = sin(pkin(6));
t158 = t144 * t145;
t146 = cos(pkin(10));
t157 = t145 * t146;
t149 = sin(qJ(2));
t156 = t145 * t149;
t147 = cos(pkin(6));
t155 = t147 * t149;
t151 = cos(qJ(2));
t154 = t147 * t151;
t153 = t148 * t151;
t152 = t150 * t151;
t141 = sin(t143);
t139 = -t144 * t155 + t146 * t151;
t138 = t144 * t154 + t146 * t149;
t137 = t144 * t151 + t146 * t155;
t136 = t144 * t149 - t146 * t154;
t135 = t147 * t141 + t142 * t156;
t134 = -t141 * t156 + t147 * t142;
t133 = t139 * t142 + t141 * t158;
t132 = -t139 * t141 + t142 * t158;
t131 = t137 * t142 - t141 * t157;
t130 = -t137 * t141 - t142 * t157;
t1 = [0, -t138 * t159 + t139 * t148, 0, t132 * t150, -t133 * t148 + t138 * t150, 0; 0, -t136 * t159 + t137 * t148, 0, t130 * t150, -t131 * t148 + t136 * t150, 0; 0 (t142 * t152 + t148 * t149) * t145, 0, t134 * t150, -t135 * t148 - t145 * t152, 0; 0, -t138 * t141, 0, t133, 0, 0; 0, -t136 * t141, 0, t131, 0, 0; 0, t145 * t151 * t141, 0, t135, 0, 0; 0, -t138 * t160 - t139 * t150, 0, t132 * t148, t133 * t150 + t138 * t148, 0; 0, -t136 * t160 - t137 * t150, 0, t130 * t148, t131 * t150 + t136 * t148, 0; 0 (t142 * t153 - t149 * t150) * t145, 0, t134 * t148, t135 * t150 - t145 * t153, 0;];
JR_rot  = t1;
