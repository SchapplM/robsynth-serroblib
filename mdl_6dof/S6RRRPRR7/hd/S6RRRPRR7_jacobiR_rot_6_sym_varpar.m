% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:07
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:07:04
% EndTime: 2019-02-22 12:07:04
% DurationCPUTime: 0.15s
% Computational Cost: add. (230->35), mult. (270->63), div. (0->0), fcn. (395->10), ass. (0->43)
t168 = cos(pkin(6));
t170 = sin(qJ(2));
t174 = cos(qJ(1));
t176 = t174 * t170;
t171 = sin(qJ(1));
t173 = cos(qJ(2));
t178 = t171 * t173;
t159 = t168 * t176 + t178;
t166 = qJ(3) + pkin(12) + qJ(5);
t164 = sin(t166);
t165 = cos(t166);
t167 = sin(pkin(6));
t181 = t167 * t174;
t152 = -t159 * t165 + t164 * t181;
t175 = t174 * t173;
t179 = t171 * t170;
t158 = -t168 * t175 + t179;
t169 = sin(qJ(6));
t172 = cos(qJ(6));
t192 = t152 * t169 + t158 * t172;
t191 = t152 * t172 - t158 * t169;
t150 = -t159 * t164 - t165 * t181;
t190 = t150 * t169;
t161 = -t168 * t179 + t175;
t182 = t167 * t171;
t153 = t161 * t164 - t165 * t182;
t189 = t153 * t169;
t183 = t167 * t170;
t156 = -t164 * t183 + t168 * t165;
t188 = t156 * t169;
t185 = t165 * t169;
t184 = t165 * t172;
t180 = t169 * t173;
t177 = t172 * t173;
t160 = t168 * t178 + t176;
t157 = t168 * t164 + t165 * t183;
t155 = t156 * t172;
t154 = t161 * t165 + t164 * t182;
t149 = t153 * t172;
t148 = t150 * t172;
t147 = t154 * t172 + t160 * t169;
t146 = -t154 * t169 + t160 * t172;
t1 = [t191, -t160 * t184 + t161 * t169, -t149, 0, -t149, t146; t147, -t158 * t184 + t159 * t169, t148, 0, t148, t192; 0 (t165 * t177 + t169 * t170) * t167, t155, 0, t155, -t157 * t169 - t167 * t177; -t192, t160 * t185 + t161 * t172, t189, 0, t189, -t147; t146, t158 * t185 + t159 * t172, -t190, 0, -t190, t191; 0 (-t165 * t180 + t170 * t172) * t167, -t188, 0, -t188, -t157 * t172 + t167 * t180; t150, -t160 * t164, t154, 0, t154, 0; t153, -t158 * t164, -t152, 0, -t152, 0; 0, t167 * t173 * t164, t157, 0, t157, 0;];
JR_rot  = t1;
