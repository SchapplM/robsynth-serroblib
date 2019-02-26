% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:37
% EndTime: 2019-02-26 21:40:37
% DurationCPUTime: 0.12s
% Computational Cost: add. (151->32), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->37)
t167 = sin(qJ(4));
t171 = cos(qJ(4));
t165 = cos(pkin(6));
t162 = sin(pkin(11));
t164 = cos(pkin(11));
t168 = sin(qJ(2));
t172 = cos(qJ(2));
t176 = t172 * t162 + t168 * t164;
t157 = t176 * t165;
t158 = t168 * t162 - t172 * t164;
t169 = sin(qJ(1));
t173 = cos(qJ(1));
t178 = t173 * t157 - t169 * t158;
t163 = sin(pkin(6));
t183 = t163 * t173;
t141 = t167 * t178 + t171 * t183;
t174 = t158 * t165;
t147 = -t169 * t176 - t173 * t174;
t166 = sin(qJ(6));
t170 = cos(qJ(6));
t188 = -t141 * t166 + t147 * t170;
t187 = t141 * t170 + t147 * t166;
t184 = t163 * t169;
t182 = t166 * t167;
t181 = t167 * t170;
t177 = -t169 * t157 - t173 * t158;
t175 = t167 * t183 - t171 * t178;
t156 = t176 * t163;
t155 = t158 * t163;
t153 = t156 * t171 + t165 * t167;
t152 = t156 * t167 - t165 * t171;
t150 = t169 * t174 - t173 * t176;
t145 = t167 * t184 + t171 * t177;
t144 = t167 * t177 - t171 * t184;
t140 = t144 * t166 - t150 * t170;
t139 = t144 * t170 + t150 * t166;
t1 = [t188, t150 * t182 + t170 * t177, 0, t145 * t166, 0, t139; t140, t147 * t182 + t170 * t178, 0, -t175 * t166, 0, t187; 0, -t155 * t182 + t156 * t170, 0, t153 * t166, 0, t152 * t170 - t155 * t166; -t187, t150 * t181 - t166 * t177, 0, t145 * t170, 0, -t140; t139, t147 * t181 - t166 * t178, 0, -t175 * t170, 0, t188; 0, -t155 * t181 - t156 * t166, 0, t153 * t170, 0, -t152 * t166 - t155 * t170; t175, t150 * t171, 0, -t144, 0, 0; t145, t147 * t171, 0, -t141, 0, 0; 0, -t155 * t171, 0, -t152, 0, 0;];
JR_rot  = t1;
