% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:46
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR12_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:46:32
% EndTime: 2019-02-22 11:46:32
% DurationCPUTime: 0.12s
% Computational Cost: add. (164->36), mult. (270->63), div. (0->0), fcn. (395->10), ass. (0->43)
t164 = cos(pkin(6));
t169 = cos(qJ(2));
t170 = cos(qJ(1));
t171 = t170 * t169;
t166 = sin(qJ(2));
t167 = sin(qJ(1));
t174 = t167 * t166;
t154 = -t164 * t171 + t174;
t162 = qJ(4) + qJ(5);
t160 = sin(t162);
t161 = cos(t162);
t163 = sin(pkin(6));
t177 = t163 * t170;
t148 = -t154 * t160 + t161 * t177;
t172 = t170 * t166;
t173 = t167 * t169;
t155 = t164 * t172 + t173;
t165 = sin(qJ(6));
t168 = cos(qJ(6));
t188 = t148 * t165 + t155 * t168;
t187 = t148 * t168 - t155 * t165;
t156 = t164 * t173 + t172;
t179 = t163 * t167;
t145 = -t156 * t161 + t160 * t179;
t186 = t145 * t165;
t147 = t154 * t161 + t160 * t177;
t185 = t147 * t165;
t178 = t163 * t169;
t152 = -t164 * t160 - t161 * t178;
t184 = t152 * t165;
t181 = t160 * t165;
t180 = t160 * t168;
t176 = t165 * t166;
t175 = t166 * t168;
t157 = -t164 * t174 + t171;
t153 = -t160 * t178 + t164 * t161;
t150 = t152 * t168;
t146 = t156 * t160 + t161 * t179;
t144 = t147 * t168;
t143 = t145 * t168;
t142 = t146 * t168 + t157 * t165;
t141 = -t146 * t165 + t157 * t168;
t1 = [t187, -t156 * t165 + t157 * t180, 0, -t143, -t143, t141; t142, -t154 * t165 + t155 * t180, 0, t144, t144, t188; 0 (t160 * t175 + t165 * t169) * t163, 0, t150, t150, -t153 * t165 + t163 * t175; -t188, -t156 * t168 - t157 * t181, 0, t186, t186, -t142; t141, -t154 * t168 - t155 * t181, 0, -t185, -t185, t187; 0 (-t160 * t176 + t168 * t169) * t163, 0, -t184, -t184, -t153 * t168 - t163 * t176; t147, -t157 * t161, 0, t146, t146, 0; t145, -t155 * t161, 0, -t148, -t148, 0; 0, -t163 * t166 * t161, 0, t153, t153, 0;];
JR_rot  = t1;
