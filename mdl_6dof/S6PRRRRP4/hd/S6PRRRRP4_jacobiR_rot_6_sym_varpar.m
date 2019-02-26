% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRP4_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:54
% EndTime: 2019-02-26 20:16:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (113->25), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->36)
t168 = qJ(4) + qJ(5);
t166 = sin(t168);
t175 = cos(qJ(3));
t184 = t166 * t175;
t167 = cos(t168);
t183 = t167 * t175;
t170 = sin(pkin(6));
t173 = sin(qJ(3));
t182 = t170 * t173;
t181 = t170 * t175;
t176 = cos(qJ(2));
t180 = t170 * t176;
t172 = cos(pkin(6));
t174 = sin(qJ(2));
t179 = t172 * t174;
t178 = t172 * t176;
t177 = t175 * t176;
t171 = cos(pkin(11));
t169 = sin(pkin(11));
t164 = t172 * t173 + t174 * t181;
t163 = t172 * t175 - t174 * t182;
t162 = -t169 * t179 + t171 * t176;
t161 = t169 * t178 + t171 * t174;
t160 = t169 * t176 + t171 * t179;
t159 = t169 * t174 - t171 * t178;
t158 = t162 * t175 + t169 * t182;
t157 = -t162 * t173 + t169 * t181;
t156 = t160 * t175 - t171 * t182;
t155 = -t160 * t173 - t171 * t181;
t154 = t164 * t167 - t166 * t180;
t153 = -t164 * t166 - t167 * t180;
t152 = t158 * t167 + t161 * t166;
t151 = -t158 * t166 + t161 * t167;
t150 = t156 * t167 + t159 * t166;
t149 = -t156 * t166 + t159 * t167;
t1 = [0, -t161 * t183 + t162 * t166, t157 * t167, t151, t151, 0; 0, -t159 * t183 + t160 * t166, t155 * t167, t149, t149, 0; 0 (t166 * t174 + t167 * t177) * t170, t163 * t167, t153, t153, 0; 0, -t161 * t173, t158, 0, 0, 0; 0, -t159 * t173, t156, 0, 0, 0; 0, t173 * t180, t164, 0, 0, 0; 0, -t161 * t184 - t162 * t167, t157 * t166, t152, t152, 0; 0, -t159 * t184 - t160 * t167, t155 * t166, t150, t150, 0; 0 (t166 * t177 - t167 * t174) * t170, t163 * t166, t154, t154, 0;];
JR_rot  = t1;
