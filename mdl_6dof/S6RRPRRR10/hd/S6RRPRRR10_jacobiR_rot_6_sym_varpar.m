% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR10
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
% Datum: 2019-02-22 11:45
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:45:15
% EndTime: 2019-02-22 11:45:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (207->33), mult. (275->62), div. (0->0), fcn. (402->10), ass. (0->40)
t165 = cos(pkin(6));
t166 = sin(qJ(2));
t169 = cos(qJ(1));
t171 = t169 * t166;
t167 = sin(qJ(1));
t168 = cos(qJ(2));
t172 = t167 * t168;
t153 = t165 * t171 + t172;
t162 = pkin(12) + qJ(4);
t158 = sin(t162);
t159 = cos(t162);
t164 = sin(pkin(6));
t174 = t164 * t169;
t147 = -t153 * t159 + t158 * t174;
t170 = t169 * t168;
t173 = t167 * t166;
t152 = -t165 * t170 + t173;
t163 = qJ(5) + qJ(6);
t160 = sin(t163);
t161 = cos(t163);
t139 = t147 * t160 + t152 * t161;
t140 = t147 * t161 - t152 * t160;
t180 = t159 * t160;
t179 = t159 * t161;
t178 = t159 * t168;
t177 = t164 * t166;
t176 = t164 * t167;
t175 = t164 * t168;
t145 = -t153 * t158 - t159 * t174;
t155 = -t165 * t173 + t170;
t154 = t165 * t172 + t171;
t151 = t158 * t165 + t159 * t177;
t150 = -t158 * t177 + t159 * t165;
t149 = t155 * t159 + t158 * t176;
t148 = t155 * t158 - t159 * t176;
t144 = -t151 * t161 + t160 * t175;
t143 = -t151 * t160 - t161 * t175;
t142 = t149 * t161 + t154 * t160;
t141 = -t149 * t160 + t154 * t161;
t1 = [t140, -t154 * t179 + t155 * t160, 0, -t148 * t161, t141, t141; t142, -t152 * t179 + t153 * t160, 0, t145 * t161, t139, t139; 0 (t160 * t166 + t161 * t178) * t164, 0, t150 * t161, t143, t143; -t139, t154 * t180 + t155 * t161, 0, t148 * t160, -t142, -t142; t141, t152 * t180 + t153 * t161, 0, -t145 * t160, t140, t140; 0 (-t160 * t178 + t161 * t166) * t164, 0, -t150 * t160, t144, t144; t145, -t154 * t158, 0, t149, 0, 0; t148, -t152 * t158, 0, -t147, 0, 0; 0, t158 * t175, 0, t151, 0, 0;];
JR_rot  = t1;
