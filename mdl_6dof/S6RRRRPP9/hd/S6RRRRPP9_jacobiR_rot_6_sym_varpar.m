% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPP9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:16
% EndTime: 2019-02-26 22:30:16
% DurationCPUTime: 0.12s
% Computational Cost: add. (71->29), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->36)
t158 = cos(pkin(6));
t161 = sin(qJ(2));
t166 = cos(qJ(1));
t168 = t166 * t161;
t162 = sin(qJ(1));
t165 = cos(qJ(2));
t171 = t162 * t165;
t152 = t158 * t168 + t171;
t160 = sin(qJ(3));
t164 = cos(qJ(3));
t157 = sin(pkin(6));
t174 = t157 * t166;
t146 = -t152 * t164 + t160 * t174;
t167 = t166 * t165;
t172 = t162 * t161;
t151 = -t158 * t167 + t172;
t159 = sin(qJ(4));
t163 = cos(qJ(4));
t181 = t146 * t159 + t151 * t163;
t180 = t146 * t163 - t151 * t159;
t177 = t157 * t160;
t176 = t157 * t164;
t175 = t157 * t165;
t173 = t159 * t164;
t170 = t163 * t164;
t169 = t164 * t165;
t144 = -t152 * t160 - t164 * t174;
t154 = -t158 * t172 + t167;
t153 = t158 * t171 + t168;
t150 = t158 * t160 + t161 * t176;
t149 = t158 * t164 - t161 * t177;
t148 = t154 * t164 + t162 * t177;
t147 = t154 * t160 - t162 * t176;
t143 = t148 * t163 + t153 * t159;
t142 = t148 * t159 - t153 * t163;
t1 = [t144, -t153 * t160, t148, 0, 0, 0; t147, -t151 * t160, -t146, 0, 0, 0; 0, t160 * t175, t150, 0, 0, 0; t181, -t153 * t173 - t154 * t163, -t147 * t159, t143, 0, 0; t142, -t151 * t173 - t152 * t163, t144 * t159, -t180, 0, 0; 0 (t159 * t169 - t161 * t163) * t157, t149 * t159, t150 * t163 - t159 * t175, 0, 0; t180, -t153 * t170 + t154 * t159, -t147 * t163, -t142, 0, 0; t143, -t151 * t170 + t152 * t159, t144 * t163, t181, 0, 0; 0 (t159 * t161 + t163 * t169) * t157, t149 * t163, -t150 * t159 - t163 * t175, 0, 0;];
JR_rot  = t1;
