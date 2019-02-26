% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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
% Datum: 2019-02-26 22:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:50:42
% EndTime: 2019-02-26 22:50:42
% DurationCPUTime: 0.11s
% Computational Cost: add. (238->33), mult. (331->62), div. (0->0), fcn. (484->10), ass. (0->39)
t159 = cos(pkin(6));
t161 = sin(qJ(2));
t165 = cos(qJ(1));
t167 = t165 * t161;
t162 = sin(qJ(1));
t164 = cos(qJ(2));
t169 = t162 * t164;
t150 = t159 * t167 + t169;
t160 = sin(qJ(3));
t163 = cos(qJ(3));
t158 = sin(pkin(6));
t171 = t158 * t165;
t144 = -t150 * t163 + t160 * t171;
t166 = t165 * t164;
t170 = t162 * t161;
t149 = -t159 * t166 + t170;
t157 = qJ(4) + qJ(5) + qJ(6);
t155 = sin(t157);
t156 = cos(t157);
t136 = t144 * t155 + t149 * t156;
t137 = t144 * t156 - t149 * t155;
t176 = t155 * t163;
t175 = t156 * t163;
t174 = t158 * t160;
t173 = t158 * t163;
t172 = t158 * t164;
t168 = t163 * t164;
t142 = -t150 * t160 - t163 * t171;
t152 = -t159 * t170 + t166;
t151 = t159 * t169 + t167;
t148 = t159 * t160 + t161 * t173;
t147 = t159 * t163 - t161 * t174;
t146 = t152 * t163 + t162 * t174;
t145 = t152 * t160 - t162 * t173;
t141 = -t148 * t156 + t155 * t172;
t140 = -t148 * t155 - t156 * t172;
t139 = t146 * t156 + t151 * t155;
t138 = -t146 * t155 + t151 * t156;
t1 = [t137, -t151 * t175 + t152 * t155, -t145 * t156, t138, t138, t138; t139, -t149 * t175 + t150 * t155, t142 * t156, t136, t136, t136; 0 (t155 * t161 + t156 * t168) * t158, t147 * t156, t140, t140, t140; -t136, t151 * t176 + t152 * t156, t145 * t155, -t139, -t139, -t139; t138, t149 * t176 + t150 * t156, -t142 * t155, t137, t137, t137; 0 (-t155 * t168 + t156 * t161) * t158, -t147 * t155, t141, t141, t141; t142, -t151 * t160, t146, 0, 0, 0; t145, -t149 * t160, -t144, 0, 0, 0; 0, t160 * t172, t148, 0, 0, 0;];
JR_rot  = t1;
