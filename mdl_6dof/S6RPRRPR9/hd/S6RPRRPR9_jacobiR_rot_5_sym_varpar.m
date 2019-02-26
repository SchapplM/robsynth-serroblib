% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR9_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:29
% EndTime: 2019-02-26 21:05:29
% DurationCPUTime: 0.11s
% Computational Cost: add. (131->31), mult. (307->59), div. (0->0), fcn. (430->12), ass. (0->41)
t173 = cos(pkin(6));
t168 = sin(pkin(12));
t177 = cos(qJ(1));
t181 = t177 * t168;
t171 = cos(pkin(12));
t175 = sin(qJ(1));
t182 = t175 * t171;
t160 = t173 * t181 + t182;
t174 = sin(qJ(3));
t176 = cos(qJ(3));
t180 = t177 * t171;
t183 = t175 * t168;
t159 = -t173 * t180 + t183;
t169 = sin(pkin(7));
t172 = cos(pkin(7));
t170 = sin(pkin(6));
t185 = t170 * t177;
t178 = t159 * t172 + t169 * t185;
t149 = -t160 * t176 + t178 * t174;
t154 = -t159 * t169 + t172 * t185;
t167 = qJ(4) + pkin(13);
t165 = sin(t167);
t166 = cos(t167);
t192 = t149 * t165 - t154 * t166;
t191 = t149 * t166 + t154 * t165;
t187 = t169 * t173;
t186 = t170 * t175;
t184 = t172 * t176;
t179 = t169 * t186;
t147 = -t160 * t174 - t178 * t176;
t162 = -t173 * t183 + t180;
t161 = -t173 * t182 - t181;
t158 = -t170 * t171 * t169 + t173 * t172;
t156 = -t161 * t169 + t172 * t186;
t153 = t174 * t187 + (t171 * t172 * t174 + t168 * t176) * t170;
t152 = t176 * t187 + (-t168 * t174 + t171 * t184) * t170;
t151 = t162 * t176 + (t161 * t172 + t179) * t174;
t150 = -t161 * t184 + t162 * t174 - t176 * t179;
t146 = t151 * t166 + t156 * t165;
t145 = -t151 * t165 + t156 * t166;
t1 = [t191, 0, -t150 * t166, t145, 0, 0; t146, 0, t147 * t166, t192, 0, 0; 0, 0, t152 * t166, -t153 * t165 + t158 * t166, 0, 0; -t192, 0, t150 * t165, -t146, 0, 0; t145, 0, -t147 * t165, t191, 0, 0; 0, 0, -t152 * t165, -t153 * t166 - t158 * t165, 0, 0; t147, 0, t151, 0, 0, 0; t150, 0, -t149, 0, 0, 0; 0, 0, t153, 0, 0, 0;];
JR_rot  = t1;
