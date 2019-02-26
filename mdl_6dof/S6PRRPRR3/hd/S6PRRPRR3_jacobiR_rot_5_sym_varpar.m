% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR3_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:24
% EndTime: 2019-02-26 20:05:24
% DurationCPUTime: 0.15s
% Computational Cost: add. (163->42), mult. (469->95), div. (0->0), fcn. (653->14), ass. (0->46)
t197 = cos(qJ(3));
t196 = cos(pkin(13));
t174 = sin(pkin(12));
t176 = sin(pkin(6));
t195 = t174 * t176;
t175 = sin(pkin(7));
t194 = t175 * t176;
t180 = sin(qJ(5));
t193 = t175 * t180;
t183 = cos(qJ(5));
t192 = t175 * t183;
t177 = cos(pkin(12));
t191 = t176 * t177;
t178 = cos(pkin(7));
t190 = t176 * t178;
t179 = cos(pkin(6));
t182 = sin(qJ(2));
t189 = t179 * t182;
t184 = cos(qJ(2));
t188 = t179 * t184;
t187 = t182 * t194;
t173 = sin(pkin(13));
t181 = sin(qJ(3));
t186 = t197 * t173 + t181 * t196;
t185 = -t181 * t173 + t197 * t196;
t161 = t186 * t175;
t163 = t186 * t178;
t165 = -t174 * t182 + t177 * t188;
t166 = t174 * t184 + t177 * t189;
t150 = -t161 * t191 + t163 * t165 + t166 * t185;
t167 = -t174 * t188 - t177 * t182;
t168 = -t174 * t189 + t177 * t184;
t152 = t161 * t195 + t163 * t167 + t168 * t185;
t154 = t179 * t161 + (t163 * t184 + t182 * t185) * t176;
t164 = t179 * t178 - t184 * t194;
t162 = t185 * t178;
t160 = t185 * t175;
t159 = -t167 * t175 + t174 * t190;
t158 = -t165 * t175 - t177 * t190;
t157 = (-t163 * t182 + t184 * t185) * t176;
t156 = -t163 * t168 + t167 * t185;
t155 = -t163 * t166 + t165 * t185;
t153 = t179 * t160 + (t162 * t184 - t182 * t186) * t176;
t151 = t160 * t195 + t162 * t167 - t168 * t186;
t149 = -t160 * t191 + t162 * t165 - t166 * t186;
t1 = [0, t156 * t183 + t168 * t193, t151 * t183, 0, -t152 * t180 + t159 * t183, 0; 0, t155 * t183 + t166 * t193, t149 * t183, 0, -t150 * t180 + t158 * t183, 0; 0, t157 * t183 + t180 * t187, t153 * t183, 0, -t154 * t180 + t164 * t183, 0; 0, -t156 * t180 + t168 * t192, -t151 * t180, 0, -t152 * t183 - t159 * t180, 0; 0, -t155 * t180 + t166 * t192, -t149 * t180, 0, -t150 * t183 - t158 * t180, 0; 0, -t157 * t180 + t183 * t187, -t153 * t180, 0, -t154 * t183 - t164 * t180, 0; 0, t162 * t168 + t167 * t186, t152, 0, 0, 0; 0, t162 * t166 + t165 * t186, t150, 0, 0, 0; 0 (t162 * t182 + t184 * t186) * t176, t154, 0, 0, 0;];
JR_rot  = t1;
