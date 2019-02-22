% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR8
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

function JR_rot = S6RRRPRR8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:07:39
% EndTime: 2019-02-22 12:07:39
% DurationCPUTime: 0.11s
% Computational Cost: add. (207->33), mult. (275->62), div. (0->0), fcn. (402->10), ass. (0->40)
t166 = cos(pkin(6));
t167 = sin(qJ(2));
t170 = cos(qJ(1));
t172 = t170 * t167;
t168 = sin(qJ(1));
t169 = cos(qJ(2));
t173 = t168 * t169;
t154 = t166 * t172 + t173;
t163 = qJ(3) + pkin(12);
t159 = sin(t163);
t160 = cos(t163);
t165 = sin(pkin(6));
t175 = t165 * t170;
t148 = -t154 * t160 + t159 * t175;
t171 = t170 * t169;
t174 = t168 * t167;
t153 = -t166 * t171 + t174;
t164 = qJ(5) + qJ(6);
t161 = sin(t164);
t162 = cos(t164);
t140 = t148 * t161 + t153 * t162;
t141 = t148 * t162 - t153 * t161;
t181 = t160 * t161;
t180 = t160 * t162;
t179 = t160 * t169;
t178 = t165 * t167;
t177 = t165 * t168;
t176 = t165 * t169;
t146 = -t154 * t159 - t160 * t175;
t156 = -t166 * t174 + t171;
t155 = t166 * t173 + t172;
t152 = t166 * t159 + t160 * t178;
t151 = -t159 * t178 + t166 * t160;
t150 = t156 * t160 + t159 * t177;
t149 = t156 * t159 - t160 * t177;
t145 = -t152 * t162 + t161 * t176;
t144 = -t152 * t161 - t162 * t176;
t143 = t150 * t162 + t155 * t161;
t142 = -t150 * t161 + t155 * t162;
t1 = [t141, -t155 * t180 + t156 * t161, -t149 * t162, 0, t142, t142; t143, -t153 * t180 + t154 * t161, t146 * t162, 0, t140, t140; 0 (t161 * t167 + t162 * t179) * t165, t151 * t162, 0, t144, t144; -t140, t155 * t181 + t156 * t162, t149 * t161, 0, -t143, -t143; t142, t153 * t181 + t154 * t162, -t146 * t161, 0, t141, t141; 0 (-t161 * t179 + t162 * t167) * t165, -t151 * t161, 0, t145, t145; t146, -t155 * t159, t150, 0, 0, 0; t149, -t153 * t159, -t148, 0, 0, 0; 0, t159 * t176, t152, 0, 0, 0;];
JR_rot  = t1;
