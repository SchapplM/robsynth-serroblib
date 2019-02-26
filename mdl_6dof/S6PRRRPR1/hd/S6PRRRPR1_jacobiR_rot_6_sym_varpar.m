% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR1
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
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:39
% EndTime: 2019-02-26 20:10:39
% DurationCPUTime: 0.12s
% Computational Cost: add. (186->31), mult. (214->65), div. (0->0), fcn. (313->10), ass. (0->37)
t149 = sin(pkin(11));
t151 = cos(pkin(11));
t156 = cos(qJ(2));
t152 = cos(pkin(6));
t154 = sin(qJ(2));
t160 = t152 * t154;
t142 = t149 * t156 + t151 * t160;
t148 = qJ(3) + qJ(4) + pkin(12);
t146 = sin(t148);
t147 = cos(t148);
t150 = sin(pkin(6));
t162 = t150 * t151;
t134 = -t142 * t146 - t147 * t162;
t153 = sin(qJ(6));
t168 = t134 * t153;
t144 = -t149 * t160 + t151 * t156;
t163 = t149 * t150;
t136 = -t144 * t146 + t147 * t163;
t167 = t136 * t153;
t161 = t150 * t154;
t139 = -t146 * t161 + t152 * t147;
t166 = t139 * t153;
t165 = t147 * t153;
t155 = cos(qJ(6));
t164 = t147 * t155;
t159 = t152 * t156;
t158 = t153 * t156;
t157 = t155 * t156;
t143 = t149 * t159 + t151 * t154;
t141 = t149 * t154 - t151 * t159;
t140 = t152 * t146 + t147 * t161;
t138 = t139 * t155;
t137 = t144 * t147 + t146 * t163;
t135 = t142 * t147 - t146 * t162;
t133 = t136 * t155;
t132 = t134 * t155;
t1 = [0, -t143 * t164 + t144 * t153, t133, t133, 0, -t137 * t153 + t143 * t155; 0, -t141 * t164 + t142 * t153, t132, t132, 0, -t135 * t153 + t141 * t155; 0 (t147 * t157 + t153 * t154) * t150, t138, t138, 0, -t140 * t153 - t150 * t157; 0, t143 * t165 + t144 * t155, -t167, -t167, 0, -t137 * t155 - t143 * t153; 0, t141 * t165 + t142 * t155, -t168, -t168, 0, -t135 * t155 - t141 * t153; 0 (-t147 * t158 + t154 * t155) * t150, -t166, -t166, 0, -t140 * t155 + t150 * t158; 0, -t143 * t146, t137, t137, 0, 0; 0, -t141 * t146, t135, t135, 0, 0; 0, t150 * t156 * t146, t140, t140, 0, 0;];
JR_rot  = t1;
