% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:28
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRR10_jacobiR_rot_3_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobiR_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobiR_rot_3_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:17
% EndTime: 2018-11-23 11:27:17
% DurationCPUTime: 0.12s
% Computational Cost: add. (248->40), mult. (241->59), div. (0->0), fcn. (255->18), ass. (0->41)
t160 = sin(pkin(6));
t165 = sin(qJ(1));
t186 = t160 * t165;
t168 = cos(qJ(1));
t185 = t160 * t168;
t184 = pkin(6) - qJ(2);
t183 = pkin(6) + qJ(2);
t182 = pkin(7) - qJ(3);
t181 = pkin(7) + qJ(3);
t180 = cos(t184);
t179 = cos(t182);
t178 = sin(t184);
t177 = sin(t182);
t176 = cos(t183) / 0.2e1;
t175 = cos(t181) / 0.2e1;
t174 = sin(t183) / 0.2e1;
t173 = sin(t181) / 0.2e1;
t148 = t174 - t178 / 0.2e1;
t167 = cos(qJ(2));
t139 = t168 * t148 + t165 * t167;
t143 = t165 * t148 - t168 * t167;
t164 = sin(qJ(2));
t169 = t180 / 0.2e1 + t176;
t138 = t165 * t164 - t168 * t169;
t145 = t173 + t177 / 0.2e1;
t150 = t179 / 0.2e1 + t175;
t163 = sin(qJ(3));
t172 = t138 * t150 + t139 * t163 + t145 * t185;
t146 = t173 - t177 / 0.2e1;
t149 = t175 - t179 / 0.2e1;
t166 = cos(qJ(3));
t171 = t138 * t146 - t139 * t166 - t149 * t185;
t141 = -t168 * t164 - t165 * t169;
t170 = -t141 * t146 + t143 * t166 + t149 * t186;
t162 = cos(pkin(6));
t161 = cos(pkin(7));
t159 = sin(pkin(7));
t151 = t176 - t180 / 0.2e1;
t147 = t174 + t178 / 0.2e1;
t137 = t141 * t150 + t143 * t163 + t145 * t186;
t1 = [t171, t141 * t166 + t143 * t146, t137, 0, 0, 0; -t170, -t138 * t166 - t139 * t146, -t172, 0, 0, 0; 0, t151 * t146 + t147 * t166, t162 * t145 + t147 * t150 + t151 * t163, 0, 0, 0; t172, -t141 * t163 + t143 * t150, t170, 0, 0, 0; t137, t138 * t163 - t139 * t150, t171, 0, 0, 0; 0, -t147 * t163 + t151 * t150, -t147 * t146 + t162 * t149 + t151 * t166, 0, 0, 0; -t138 * t159 + t161 * t185, -t143 * t159, 0, 0, 0, 0; -t141 * t159 + t161 * t186, t139 * t159, 0, 0, 0, 0; 0, -t151 * t159, 0, 0, 0, 0;];
JR_rot  = t1;
