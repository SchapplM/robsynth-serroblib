% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:44
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:44:42
% EndTime: 2019-02-22 11:44:42
% DurationCPUTime: 0.14s
% Computational Cost: add. (230->35), mult. (270->63), div. (0->0), fcn. (395->10), ass. (0->43)
t167 = cos(pkin(6));
t169 = sin(qJ(2));
t173 = cos(qJ(1));
t175 = t173 * t169;
t170 = sin(qJ(1));
t172 = cos(qJ(2));
t177 = t170 * t172;
t158 = t167 * t175 + t177;
t165 = pkin(12) + qJ(4) + qJ(5);
t163 = sin(t165);
t164 = cos(t165);
t166 = sin(pkin(6));
t180 = t166 * t173;
t151 = -t158 * t164 + t163 * t180;
t174 = t173 * t172;
t178 = t170 * t169;
t157 = -t167 * t174 + t178;
t168 = sin(qJ(6));
t171 = cos(qJ(6));
t191 = t151 * t168 + t157 * t171;
t190 = t151 * t171 - t157 * t168;
t149 = -t158 * t163 - t164 * t180;
t189 = t149 * t168;
t160 = -t167 * t178 + t174;
t181 = t166 * t170;
t152 = t160 * t163 - t164 * t181;
t188 = t152 * t168;
t182 = t166 * t169;
t155 = -t163 * t182 + t167 * t164;
t187 = t155 * t168;
t184 = t164 * t168;
t183 = t164 * t171;
t179 = t168 * t172;
t176 = t171 * t172;
t159 = t167 * t177 + t175;
t156 = t167 * t163 + t164 * t182;
t154 = t155 * t171;
t153 = t160 * t164 + t163 * t181;
t148 = t152 * t171;
t147 = t149 * t171;
t146 = t153 * t171 + t159 * t168;
t145 = -t153 * t168 + t159 * t171;
t1 = [t190, -t159 * t183 + t160 * t168, 0, -t148, -t148, t145; t146, -t157 * t183 + t158 * t168, 0, t147, t147, t191; 0 (t164 * t176 + t168 * t169) * t166, 0, t154, t154, -t156 * t168 - t166 * t176; -t191, t159 * t184 + t160 * t171, 0, t188, t188, -t146; t145, t157 * t184 + t158 * t171, 0, -t189, -t189, t190; 0 (-t164 * t179 + t169 * t171) * t166, 0, -t187, -t187, -t156 * t171 + t166 * t179; t149, -t159 * t163, 0, t153, t153, 0; t152, -t157 * t163, 0, -t151, -t151, 0; 0, t166 * t172 * t163, 0, t156, t156, 0;];
JR_rot  = t1;
