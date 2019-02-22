% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR13
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
% Datum: 2019-02-22 12:24
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR13_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:24:35
% EndTime: 2019-02-22 12:24:35
% DurationCPUTime: 0.14s
% Computational Cost: add. (71->29), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->36)
t164 = cos(pkin(6));
t167 = sin(qJ(2));
t172 = cos(qJ(1));
t174 = t172 * t167;
t168 = sin(qJ(1));
t171 = cos(qJ(2));
t177 = t168 * t171;
t158 = t164 * t174 + t177;
t166 = sin(qJ(3));
t170 = cos(qJ(3));
t163 = sin(pkin(6));
t180 = t163 * t172;
t152 = -t158 * t170 + t166 * t180;
t173 = t172 * t171;
t178 = t168 * t167;
t157 = -t164 * t173 + t178;
t165 = sin(qJ(4));
t169 = cos(qJ(4));
t187 = t152 * t165 + t157 * t169;
t186 = t152 * t169 - t157 * t165;
t183 = t163 * t166;
t182 = t163 * t170;
t181 = t163 * t171;
t179 = t165 * t170;
t176 = t169 * t170;
t175 = t170 * t171;
t150 = -t158 * t166 - t170 * t180;
t160 = -t164 * t178 + t173;
t159 = t164 * t177 + t174;
t156 = t164 * t166 + t167 * t182;
t155 = t164 * t170 - t167 * t183;
t154 = t160 * t170 + t168 * t183;
t153 = t160 * t166 - t168 * t182;
t149 = t154 * t169 + t159 * t165;
t148 = t154 * t165 - t159 * t169;
t1 = [t186, -t159 * t176 + t160 * t165, -t153 * t169, -t148, 0, 0; t149, -t157 * t176 + t158 * t165, t150 * t169, t187, 0, 0; 0 (t165 * t167 + t169 * t175) * t163, t155 * t169, -t156 * t165 - t169 * t181, 0, 0; t150, -t159 * t166, t154, 0, 0, 0; t153, -t157 * t166, -t152, 0, 0, 0; 0, t166 * t181, t156, 0, 0, 0; t187, -t159 * t179 - t160 * t169, -t153 * t165, t149, 0, 0; t148, -t157 * t179 - t158 * t169, t150 * t165, -t186, 0, 0; 0 (t165 * t175 - t167 * t169) * t163, t155 * t165, t156 * t169 - t165 * t181, 0, 0;];
JR_rot  = t1;
