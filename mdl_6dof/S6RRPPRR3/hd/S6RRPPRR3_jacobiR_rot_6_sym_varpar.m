% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:34
% EndTime: 2019-02-26 21:29:34
% DurationCPUTime: 0.16s
% Computational Cost: add. (205->32), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->38)
t176 = pkin(12) + qJ(5);
t174 = sin(t176);
t175 = cos(t176);
t180 = cos(pkin(6));
t177 = sin(pkin(11));
t179 = cos(pkin(11));
t182 = sin(qJ(2));
t185 = cos(qJ(2));
t188 = t177 * t185 + t179 * t182;
t168 = t188 * t180;
t169 = t177 * t182 - t185 * t179;
t183 = sin(qJ(1));
t186 = cos(qJ(1));
t190 = t168 * t186 - t183 * t169;
t178 = sin(pkin(6));
t192 = t178 * t186;
t154 = t174 * t192 - t175 * t190;
t187 = t169 * t180;
t158 = -t183 * t188 - t186 * t187;
t181 = sin(qJ(6));
t184 = cos(qJ(6));
t200 = t154 * t181 - t158 * t184;
t199 = t154 * t184 + t158 * t181;
t196 = t175 * t181;
t195 = t175 * t184;
t193 = t178 * t183;
t189 = -t183 * t168 - t169 * t186;
t152 = -t174 * t190 - t175 * t192;
t167 = t188 * t178;
t166 = t169 * t178;
t164 = t167 * t175 + t174 * t180;
t163 = -t167 * t174 + t175 * t180;
t161 = t183 * t187 - t186 * t188;
t156 = t174 * t193 + t175 * t189;
t155 = t174 * t189 - t175 * t193;
t151 = t156 * t184 - t161 * t181;
t150 = -t156 * t181 - t161 * t184;
t1 = [t199, t161 * t195 + t181 * t189, 0, 0, -t155 * t184, t150; t151, t158 * t195 + t181 * t190, 0, 0, t152 * t184, t200; 0, -t166 * t195 + t167 * t181, 0, 0, t163 * t184, -t164 * t181 + t166 * t184; -t200, -t161 * t196 + t184 * t189, 0, 0, t155 * t181, -t151; t150, -t158 * t196 + t184 * t190, 0, 0, -t152 * t181, t199; 0, t166 * t196 + t167 * t184, 0, 0, -t163 * t181, -t164 * t184 - t166 * t181; t152, t161 * t174, 0, 0, t156, 0; t155, t158 * t174, 0, 0, -t154, 0; 0, -t166 * t174, 0, 0, t164, 0;];
JR_rot  = t1;
