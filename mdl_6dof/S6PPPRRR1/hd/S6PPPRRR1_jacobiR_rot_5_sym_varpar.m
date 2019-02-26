% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPPRRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiR_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:44
% EndTime: 2019-02-26 19:38:44
% DurationCPUTime: 0.14s
% Computational Cost: add. (200->38), mult. (582->82), div. (0->0), fcn. (794->16), ass. (0->50)
t177 = sin(pkin(12));
t186 = cos(pkin(6));
t201 = t177 * t186;
t179 = sin(pkin(7));
t180 = sin(pkin(6));
t200 = t179 * t180;
t199 = t179 * t186;
t185 = cos(pkin(7));
t198 = t180 * t185;
t182 = cos(pkin(13));
t197 = t182 * t185;
t183 = cos(pkin(12));
t196 = t183 * t186;
t176 = sin(pkin(13));
t172 = t176 * t196 + t177 * t182;
t175 = sin(pkin(14));
t181 = cos(pkin(14));
t171 = -t177 * t176 + t182 * t196;
t192 = t171 * t185 - t183 * t200;
t161 = -t172 * t175 + t192 * t181;
t168 = -t171 * t179 - t183 * t198;
t178 = sin(pkin(8));
t184 = cos(pkin(8));
t195 = t161 * t184 + t168 * t178;
t174 = -t176 * t201 + t183 * t182;
t173 = -t183 * t176 - t182 * t201;
t191 = t173 * t185 + t177 * t200;
t163 = -t174 * t175 + t191 * t181;
t169 = -t173 * t179 + t177 * t198;
t194 = t163 * t184 + t169 * t178;
t166 = t181 * t199 + (-t175 * t176 + t181 * t197) * t180;
t170 = -t182 * t200 + t186 * t185;
t193 = t166 * t184 + t170 * t178;
t190 = cos(qJ(4));
t189 = cos(qJ(5));
t188 = sin(qJ(4));
t187 = sin(qJ(5));
t167 = t180 * t176 * t181 + (t180 * t197 + t199) * t175;
t165 = -t166 * t178 + t170 * t184;
t164 = t174 * t181 + t191 * t175;
t162 = t172 * t181 + t192 * t175;
t160 = -t163 * t178 + t169 * t184;
t159 = -t161 * t178 + t168 * t184;
t158 = t167 * t190 + t193 * t188;
t157 = -t167 * t188 + t193 * t190;
t156 = t164 * t190 + t194 * t188;
t155 = -t164 * t188 + t194 * t190;
t154 = t162 * t190 + t195 * t188;
t153 = -t162 * t188 + t195 * t190;
t1 = [0, 0, 0, t155 * t189, -t156 * t187 + t160 * t189, 0; 0, 0, 0, t153 * t189, -t154 * t187 + t159 * t189, 0; 0, 0, 0, t157 * t189, -t158 * t187 + t165 * t189, 0; 0, 0, 0, -t155 * t187, -t156 * t189 - t160 * t187, 0; 0, 0, 0, -t153 * t187, -t154 * t189 - t159 * t187, 0; 0, 0, 0, -t157 * t187, -t158 * t189 - t165 * t187, 0; 0, 0, 0, t156, 0, 0; 0, 0, 0, t154, 0, 0; 0, 0, 0, t158, 0, 0;];
JR_rot  = t1;
