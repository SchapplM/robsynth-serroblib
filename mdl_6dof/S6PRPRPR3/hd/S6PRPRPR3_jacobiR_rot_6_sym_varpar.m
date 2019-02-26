% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRPR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:47:39
% EndTime: 2019-02-26 19:47:39
% DurationCPUTime: 0.09s
% Computational Cost: add. (111->28), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->33)
t147 = sin(pkin(11));
t150 = cos(pkin(11));
t155 = sin(qJ(2));
t158 = cos(qJ(2));
t143 = t155 * t147 - t158 * t150;
t149 = sin(pkin(6));
t154 = sin(qJ(4));
t168 = t149 * t154;
t157 = cos(qJ(4));
t167 = t149 * t157;
t153 = sin(qJ(6));
t166 = t153 * t154;
t156 = cos(qJ(6));
t165 = t154 * t156;
t152 = cos(pkin(6));
t160 = t158 * t147 + t155 * t150;
t142 = t160 * t152;
t148 = sin(pkin(10));
t151 = cos(pkin(10));
t162 = t151 * t142 - t148 * t143;
t161 = -t148 * t142 - t151 * t143;
t159 = t143 * t152;
t141 = t160 * t149;
t140 = t143 * t149;
t138 = t141 * t157 + t152 * t154;
t137 = t141 * t154 - t152 * t157;
t135 = t148 * t159 - t151 * t160;
t132 = -t148 * t160 - t151 * t159;
t130 = t148 * t168 + t157 * t161;
t129 = -t148 * t167 + t154 * t161;
t128 = -t151 * t168 + t157 * t162;
t127 = t151 * t167 + t154 * t162;
t1 = [0, t135 * t166 + t156 * t161, 0, t130 * t153, 0, t129 * t156 + t135 * t153; 0, t132 * t166 + t156 * t162, 0, t128 * t153, 0, t127 * t156 + t132 * t153; 0, -t140 * t166 + t141 * t156, 0, t138 * t153, 0, t137 * t156 - t140 * t153; 0, t135 * t165 - t153 * t161, 0, t130 * t156, 0, -t129 * t153 + t135 * t156; 0, t132 * t165 - t153 * t162, 0, t128 * t156, 0, -t127 * t153 + t132 * t156; 0, -t140 * t165 - t141 * t153, 0, t138 * t156, 0, -t137 * t153 - t140 * t156; 0, t135 * t157, 0, -t129, 0, 0; 0, t132 * t157, 0, -t127, 0, 0; 0, -t140 * t157, 0, -t137, 0, 0;];
JR_rot  = t1;
