% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:50
% EndTime: 2019-02-26 22:05:51
% DurationCPUTime: 0.10s
% Computational Cost: add. (163->32), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->38)
t147 = cos(pkin(6));
t148 = sin(qJ(2));
t151 = cos(qJ(1));
t153 = t151 * t148;
t149 = sin(qJ(1));
t150 = cos(qJ(2));
t154 = t149 * t150;
t135 = t147 * t153 + t154;
t145 = qJ(3) + pkin(11);
t141 = sin(t145);
t143 = cos(t145);
t146 = sin(pkin(6));
t156 = t146 * t151;
t129 = -t135 * t143 + t141 * t156;
t152 = t151 * t150;
t155 = t149 * t148;
t134 = -t147 * t152 + t155;
t144 = pkin(12) + qJ(6);
t140 = sin(t144);
t142 = cos(t144);
t166 = t129 * t140 + t134 * t142;
t165 = t129 * t142 - t134 * t140;
t162 = t140 * t143;
t161 = t142 * t143;
t160 = t143 * t150;
t159 = t146 * t148;
t158 = t146 * t149;
t157 = t146 * t150;
t127 = -t135 * t141 - t143 * t156;
t137 = -t147 * t155 + t152;
t136 = t147 * t154 + t153;
t133 = t147 * t141 + t143 * t159;
t132 = -t141 * t159 + t147 * t143;
t131 = t137 * t143 + t141 * t158;
t130 = t137 * t141 - t143 * t158;
t126 = t131 * t142 + t136 * t140;
t125 = -t131 * t140 + t136 * t142;
t1 = [t165, -t136 * t161 + t137 * t140, -t130 * t142, 0, 0, t125; t126, -t134 * t161 + t135 * t140, t127 * t142, 0, 0, t166; 0 (t140 * t148 + t142 * t160) * t146, t132 * t142, 0, 0, -t133 * t140 - t142 * t157; -t166, t136 * t162 + t137 * t142, t130 * t140, 0, 0, -t126; t125, t134 * t162 + t135 * t142, -t127 * t140, 0, 0, t165; 0 (-t140 * t160 + t142 * t148) * t146, -t132 * t140, 0, 0, -t133 * t142 + t140 * t157; t127, -t136 * t141, t131, 0, 0, 0; t130, -t134 * t141, -t129, 0, 0, 0; 0, t141 * t157, t133, 0, 0, 0;];
JR_rot  = t1;
