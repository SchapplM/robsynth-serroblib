% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:03
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR4_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:03:49
% EndTime: 2019-02-22 10:03:49
% DurationCPUTime: 0.17s
% Computational Cost: add. (178->39), mult. (408->88), div. (0->0), fcn. (571->12), ass. (0->52)
t176 = qJ(4) + qJ(5);
t174 = sin(t176);
t178 = sin(pkin(7));
t202 = t174 * t178;
t175 = cos(t176);
t201 = t175 * t178;
t179 = sin(pkin(6));
t200 = t178 * t179;
t182 = cos(pkin(6));
t199 = t178 * t182;
t181 = cos(pkin(7));
t198 = t179 * t181;
t183 = sin(qJ(3));
t197 = t181 * t183;
t185 = cos(qJ(3));
t196 = t181 * t185;
t184 = sin(qJ(2));
t195 = t182 * t184;
t186 = cos(qJ(2));
t194 = t182 * t186;
t193 = t183 * t184;
t192 = t183 * t186;
t191 = t184 * t185;
t190 = t185 * t186;
t189 = t184 * t200;
t177 = sin(pkin(13));
t180 = cos(pkin(13));
t169 = -t177 * t184 + t180 * t194;
t188 = t169 * t181 - t180 * t200;
t171 = -t177 * t194 - t180 * t184;
t187 = t171 * t181 + t177 * t200;
t172 = -t177 * t195 + t180 * t186;
t170 = t177 * t186 + t180 * t195;
t168 = t182 * t181 - t186 * t200;
t167 = (-t181 * t193 + t190) * t179;
t166 = -t171 * t178 + t177 * t198;
t165 = -t169 * t178 - t180 * t198;
t164 = t183 * t199 + (t181 * t192 + t191) * t179;
t163 = t185 * t199 + (t181 * t190 - t193) * t179;
t162 = t171 * t185 - t172 * t197;
t161 = t169 * t185 - t170 * t197;
t160 = t172 * t185 + t187 * t183;
t159 = -t172 * t183 + t187 * t185;
t158 = t170 * t185 + t188 * t183;
t157 = -t170 * t183 + t188 * t185;
t156 = -t164 * t175 - t168 * t174;
t155 = -t164 * t174 + t168 * t175;
t154 = -t160 * t175 - t166 * t174;
t153 = -t160 * t174 + t166 * t175;
t152 = -t158 * t175 - t165 * t174;
t151 = -t158 * t174 + t165 * t175;
t1 = [0, t162 * t175 + t172 * t202, t159 * t175, t153, t153, 0; 0, t161 * t175 + t170 * t202, t157 * t175, t151, t151, 0; 0, t167 * t175 + t174 * t189, t163 * t175, t155, t155, 0; 0, -t162 * t174 + t172 * t201, -t159 * t174, t154, t154, 0; 0, -t161 * t174 + t170 * t201, -t157 * t174, t152, t152, 0; 0, -t167 * t174 + t175 * t189, -t163 * t174, t156, t156, 0; 0, t171 * t183 + t172 * t196, t160, 0, 0, 0; 0, t169 * t183 + t170 * t196, t158, 0, 0, 0; 0 (t181 * t191 + t192) * t179, t164, 0, 0, 0;];
JR_rot  = t1;
