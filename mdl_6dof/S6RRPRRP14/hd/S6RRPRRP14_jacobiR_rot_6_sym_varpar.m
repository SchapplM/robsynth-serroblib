% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:39
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP14_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:39:21
% EndTime: 2019-02-22 11:39:21
% DurationCPUTime: 0.14s
% Computational Cost: add. (74->30), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->36)
t168 = cos(pkin(6));
t175 = cos(qJ(2));
t176 = cos(qJ(1));
t177 = t176 * t175;
t171 = sin(qJ(2));
t172 = sin(qJ(1));
t180 = t172 * t171;
t161 = -t168 * t177 + t180;
t170 = sin(qJ(4));
t174 = cos(qJ(4));
t167 = sin(pkin(6));
t185 = t167 * t176;
t156 = -t161 * t170 + t174 * t185;
t178 = t176 * t171;
t179 = t172 * t175;
t162 = t168 * t178 + t179;
t169 = sin(qJ(5));
t173 = cos(qJ(5));
t191 = t156 * t169 + t162 * t173;
t190 = t156 * t173 - t162 * t169;
t187 = t167 * t174;
t186 = t167 * t175;
t184 = t169 * t170;
t183 = t169 * t171;
t182 = t170 * t173;
t181 = t171 * t173;
t155 = t161 * t174 + t170 * t185;
t164 = -t168 * t180 + t177;
t163 = t168 * t179 + t178;
t160 = t168 * t174 - t170 * t186;
t159 = -t168 * t170 - t174 * t186;
t154 = t163 * t170 + t172 * t187;
t153 = t172 * t167 * t170 - t163 * t174;
t152 = t154 * t173 + t164 * t169;
t151 = t154 * t169 - t164 * t173;
t1 = [t190, -t163 * t169 + t164 * t182, 0, -t153 * t173, -t151, 0; t152, -t161 * t169 + t162 * t182, 0, t155 * t173, t191, 0; 0 (t169 * t175 + t170 * t181) * t167, 0, t159 * t173, -t160 * t169 + t167 * t181, 0; t155, -t164 * t174, 0, t154, 0, 0; t153, -t162 * t174, 0, -t156, 0, 0; 0, -t171 * t187, 0, t160, 0, 0; t191, t163 * t173 + t164 * t184, 0, -t153 * t169, t152, 0; t151, t161 * t173 + t162 * t184, 0, t155 * t169, -t190, 0; 0 (t170 * t183 - t173 * t175) * t167, 0, t159 * t169, t160 * t173 + t167 * t183, 0;];
JR_rot  = t1;
