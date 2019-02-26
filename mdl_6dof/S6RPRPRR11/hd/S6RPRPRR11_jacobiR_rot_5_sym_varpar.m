% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR11_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobiR_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:31
% EndTime: 2019-02-26 20:54:31
% DurationCPUTime: 0.11s
% Computational Cost: add. (131->31), mult. (307->59), div. (0->0), fcn. (430->12), ass. (0->41)
t170 = cos(pkin(6));
t165 = sin(pkin(12));
t174 = cos(qJ(1));
t178 = t174 * t165;
t168 = cos(pkin(12));
t172 = sin(qJ(1));
t179 = t172 * t168;
t157 = t170 * t178 + t179;
t171 = sin(qJ(3));
t173 = cos(qJ(3));
t177 = t174 * t168;
t180 = t172 * t165;
t156 = -t170 * t177 + t180;
t166 = sin(pkin(7));
t169 = cos(pkin(7));
t167 = sin(pkin(6));
t182 = t167 * t174;
t175 = t156 * t169 + t166 * t182;
t146 = -t157 * t173 + t175 * t171;
t151 = -t156 * t166 + t169 * t182;
t164 = pkin(13) + qJ(5);
t162 = sin(t164);
t163 = cos(t164);
t189 = t146 * t162 - t151 * t163;
t188 = t146 * t163 + t151 * t162;
t184 = t166 * t170;
t183 = t167 * t172;
t181 = t169 * t173;
t176 = t166 * t183;
t144 = -t157 * t171 - t175 * t173;
t159 = -t170 * t180 + t177;
t158 = -t170 * t179 - t178;
t155 = -t167 * t168 * t166 + t170 * t169;
t153 = -t158 * t166 + t169 * t183;
t150 = t171 * t184 + (t168 * t169 * t171 + t165 * t173) * t167;
t149 = t173 * t184 + (-t165 * t171 + t168 * t181) * t167;
t148 = t159 * t173 + (t158 * t169 + t176) * t171;
t147 = -t158 * t181 + t159 * t171 - t173 * t176;
t143 = t148 * t163 + t153 * t162;
t142 = -t148 * t162 + t153 * t163;
t1 = [t188, 0, -t147 * t163, 0, t142, 0; t143, 0, t144 * t163, 0, t189, 0; 0, 0, t149 * t163, 0, -t150 * t162 + t155 * t163, 0; -t189, 0, t147 * t162, 0, -t143, 0; t142, 0, -t144 * t162, 0, t188, 0; 0, 0, -t149 * t162, 0, -t150 * t163 - t155 * t162, 0; t144, 0, t148, 0, 0, 0; t147, 0, -t146, 0, 0, 0; 0, 0, t150, 0, 0, 0;];
JR_rot  = t1;
