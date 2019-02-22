% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:22
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:22:37
% EndTime: 2019-02-22 12:22:37
% DurationCPUTime: 0.16s
% Computational Cost: add. (155->37), mult. (270->65), div. (0->0), fcn. (395->10), ass. (0->41)
t169 = sin(qJ(2));
t170 = sin(qJ(1));
t172 = cos(qJ(2));
t173 = cos(qJ(1));
t186 = cos(pkin(6));
t175 = t173 * t186;
t158 = t169 * t175 + t170 * t172;
t166 = qJ(3) + qJ(4);
t164 = sin(t166);
t165 = cos(t166);
t167 = sin(pkin(6));
t179 = t167 * t173;
t148 = t158 * t164 + t165 * t179;
t157 = t170 * t169 - t172 * t175;
t168 = sin(qJ(6));
t171 = cos(qJ(6));
t188 = -t148 * t168 - t157 * t171;
t187 = t148 * t171 - t157 * t168;
t183 = t164 * t168;
t182 = t164 * t171;
t181 = t167 * t169;
t180 = t167 * t170;
t178 = t168 * t172;
t177 = t171 * t172;
t176 = t170 * t186;
t174 = -t158 * t165 + t164 * t179;
t160 = -t169 * t176 + t173 * t172;
t159 = t173 * t169 + t172 * t176;
t156 = t186 * t164 + t165 * t181;
t155 = t164 * t181 - t186 * t165;
t154 = t156 * t171;
t153 = t156 * t168;
t152 = t160 * t165 + t164 * t180;
t151 = t160 * t164 - t165 * t180;
t147 = t152 * t171;
t146 = t152 * t168;
t145 = t174 * t171;
t144 = t174 * t168;
t143 = t151 * t168 + t159 * t171;
t142 = t151 * t171 - t159 * t168;
t1 = [t188, -t159 * t183 + t160 * t171, t146, t146, 0, t142; t143, -t157 * t183 + t158 * t171, -t144, -t144, 0, t187; 0 (t164 * t178 + t169 * t171) * t167, t153, t153, 0, t155 * t171 + t167 * t178; -t187, -t159 * t182 - t160 * t168, t147, t147, 0, -t143; t142, -t157 * t182 - t158 * t168, -t145, -t145, 0, t188; 0 (t164 * t177 - t168 * t169) * t167, t154, t154, 0, -t155 * t168 + t167 * t177; t174, -t159 * t165, -t151, -t151, 0, 0; t152, -t157 * t165, -t148, -t148, 0, 0; 0, t167 * t172 * t165, -t155, -t155, 0, 0;];
JR_rot  = t1;
