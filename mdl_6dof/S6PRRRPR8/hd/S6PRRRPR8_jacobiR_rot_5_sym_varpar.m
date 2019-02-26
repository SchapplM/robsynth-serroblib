% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR8_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:42
% EndTime: 2019-02-26 20:14:42
% DurationCPUTime: 0.09s
% Computational Cost: add. (100->38), mult. (304->88), div. (0->0), fcn. (425->12), ass. (0->45)
t173 = sin(pkin(7));
t174 = sin(pkin(6));
t199 = t173 * t174;
t177 = cos(pkin(6));
t198 = t173 * t177;
t178 = sin(qJ(4));
t197 = t173 * t178;
t181 = cos(qJ(4));
t196 = t173 * t181;
t176 = cos(pkin(7));
t195 = t174 * t176;
t179 = sin(qJ(3));
t194 = t176 * t179;
t182 = cos(qJ(3));
t193 = t176 * t182;
t180 = sin(qJ(2));
t192 = t177 * t180;
t183 = cos(qJ(2));
t191 = t177 * t183;
t190 = t179 * t180;
t189 = t179 * t183;
t188 = t180 * t182;
t187 = t182 * t183;
t186 = t180 * t199;
t172 = sin(pkin(12));
t175 = cos(pkin(12));
t167 = -t172 * t180 + t175 * t191;
t185 = t167 * t176 - t175 * t199;
t169 = -t172 * t191 - t175 * t180;
t184 = t169 * t176 + t172 * t199;
t170 = -t172 * t192 + t175 * t183;
t168 = t172 * t183 + t175 * t192;
t166 = t177 * t176 - t183 * t199;
t165 = (-t176 * t190 + t187) * t174;
t164 = -t169 * t173 + t172 * t195;
t163 = -t167 * t173 - t175 * t195;
t162 = t179 * t198 + (t176 * t189 + t188) * t174;
t161 = t182 * t198 + (t176 * t187 - t190) * t174;
t160 = t169 * t182 - t170 * t194;
t159 = t167 * t182 - t168 * t194;
t158 = t170 * t182 + t184 * t179;
t157 = -t170 * t179 + t184 * t182;
t156 = t168 * t182 + t185 * t179;
t155 = -t168 * t179 + t185 * t182;
t1 = [0, t169 * t179 + t170 * t193, t158, 0, 0, 0; 0, t167 * t179 + t168 * t193, t156, 0, 0, 0; 0 (t176 * t188 + t189) * t174, t162, 0, 0, 0; 0, -t160 * t181 - t170 * t197, -t157 * t181, t158 * t178 - t164 * t181, 0, 0; 0, -t159 * t181 - t168 * t197, -t155 * t181, t156 * t178 - t163 * t181, 0, 0; 0, -t165 * t181 - t178 * t186, -t161 * t181, t162 * t178 - t166 * t181, 0, 0; 0, t160 * t178 - t170 * t196, t157 * t178, t158 * t181 + t164 * t178, 0, 0; 0, t159 * t178 - t168 * t196, t155 * t178, t156 * t181 + t163 * t178, 0, 0; 0, t165 * t178 - t181 * t186, t161 * t178, t162 * t181 + t166 * t178, 0, 0;];
JR_rot  = t1;
