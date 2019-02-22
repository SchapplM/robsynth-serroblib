% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:29
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:29:51
% EndTime: 2019-02-22 12:29:51
% DurationCPUTime: 0.14s
% Computational Cost: add. (161->35), mult. (270->63), div. (0->0), fcn. (395->10), ass. (0->43)
t170 = cos(pkin(6));
t172 = sin(qJ(2));
t176 = cos(qJ(1));
t178 = t176 * t172;
t173 = sin(qJ(1));
t175 = cos(qJ(2));
t180 = t173 * t175;
t161 = t170 * t178 + t180;
t168 = qJ(3) + qJ(4);
t166 = sin(t168);
t167 = cos(t168);
t169 = sin(pkin(6));
t183 = t169 * t176;
t154 = -t161 * t167 + t166 * t183;
t177 = t176 * t175;
t181 = t173 * t172;
t160 = -t170 * t177 + t181;
t171 = sin(qJ(5));
t174 = cos(qJ(5));
t194 = t154 * t171 + t160 * t174;
t193 = t154 * t174 - t160 * t171;
t152 = -t161 * t166 - t167 * t183;
t192 = t152 * t171;
t163 = -t170 * t181 + t177;
t184 = t169 * t173;
t155 = t163 * t166 - t167 * t184;
t191 = t155 * t171;
t185 = t169 * t172;
t158 = -t166 * t185 + t170 * t167;
t190 = t158 * t171;
t187 = t167 * t171;
t186 = t167 * t174;
t182 = t171 * t175;
t179 = t174 * t175;
t162 = t170 * t180 + t178;
t159 = t170 * t166 + t167 * t185;
t157 = t158 * t174;
t156 = t163 * t167 + t166 * t184;
t151 = t155 * t174;
t150 = t152 * t174;
t149 = t156 * t174 + t162 * t171;
t148 = -t156 * t171 + t162 * t174;
t1 = [t193, -t162 * t186 + t163 * t171, -t151, -t151, t148, 0; t149, -t160 * t186 + t161 * t171, t150, t150, t194, 0; 0 (t167 * t179 + t171 * t172) * t169, t157, t157, -t159 * t171 - t169 * t179, 0; -t194, t162 * t187 + t163 * t174, t191, t191, -t149, 0; t148, t160 * t187 + t161 * t174, -t192, -t192, t193, 0; 0 (-t167 * t182 + t172 * t174) * t169, -t190, -t190, -t159 * t174 + t169 * t182, 0; t152, -t162 * t166, t156, t156, 0, 0; t155, -t160 * t166, -t154, -t154, 0, 0; 0, t169 * t175 * t166, t159, t159, 0, 0;];
JR_rot  = t1;
