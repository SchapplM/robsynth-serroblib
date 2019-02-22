% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:37
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:37:08
% EndTime: 2019-02-22 12:37:08
% DurationCPUTime: 0.11s
% Computational Cost: add. (144->32), mult. (275->62), div. (0->0), fcn. (402->10), ass. (0->39)
t154 = cos(pkin(6));
t156 = sin(qJ(2));
t160 = cos(qJ(1));
t162 = t160 * t156;
t157 = sin(qJ(1));
t159 = cos(qJ(2));
t164 = t157 * t159;
t145 = t154 * t162 + t164;
t155 = sin(qJ(3));
t158 = cos(qJ(3));
t153 = sin(pkin(6));
t166 = t153 * t160;
t139 = -t145 * t158 + t155 * t166;
t161 = t160 * t159;
t165 = t157 * t156;
t144 = -t154 * t161 + t165;
t152 = qJ(4) + qJ(5);
t150 = sin(t152);
t151 = cos(t152);
t131 = t139 * t150 + t144 * t151;
t132 = t139 * t151 - t144 * t150;
t171 = t150 * t158;
t170 = t151 * t158;
t169 = t153 * t155;
t168 = t153 * t158;
t167 = t153 * t159;
t163 = t158 * t159;
t137 = -t145 * t155 - t158 * t166;
t147 = -t154 * t165 + t161;
t146 = t154 * t164 + t162;
t143 = t154 * t155 + t156 * t168;
t142 = t154 * t158 - t156 * t169;
t141 = t147 * t158 + t157 * t169;
t140 = t147 * t155 - t157 * t168;
t136 = -t143 * t151 + t150 * t167;
t135 = -t143 * t150 - t151 * t167;
t134 = t141 * t151 + t146 * t150;
t133 = -t141 * t150 + t146 * t151;
t1 = [t132, -t146 * t170 + t147 * t150, -t140 * t151, t133, t133, 0; t134, -t144 * t170 + t145 * t150, t137 * t151, t131, t131, 0; 0 (t150 * t156 + t151 * t163) * t153, t142 * t151, t135, t135, 0; -t131, t146 * t171 + t147 * t151, t140 * t150, -t134, -t134, 0; t133, t144 * t171 + t145 * t151, -t137 * t150, t132, t132, 0; 0 (-t150 * t163 + t151 * t156) * t153, -t142 * t150, t136, t136, 0; t137, -t146 * t155, t141, 0, 0, 0; t140, -t144 * t155, -t139, 0, 0, 0; 0, t155 * t167, t143, 0, 0, 0;];
JR_rot  = t1;
