% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP12_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:15:32
% EndTime: 2019-02-26 22:15:32
% DurationCPUTime: 0.12s
% Computational Cost: add. (74->32), mult. (219->62), div. (0->0), fcn. (320->10), ass. (0->36)
t163 = cos(pkin(6));
t166 = sin(qJ(2));
t171 = cos(qJ(1));
t174 = t171 * t166;
t167 = sin(qJ(1));
t170 = cos(qJ(2));
t176 = t167 * t170;
t158 = t163 * t174 + t176;
t165 = sin(qJ(3));
t169 = cos(qJ(3));
t162 = sin(pkin(6));
t181 = t162 * t171;
t150 = t158 * t165 + t169 * t181;
t173 = t171 * t170;
t177 = t167 * t166;
t157 = -t163 * t173 + t177;
t164 = sin(qJ(5));
t168 = cos(qJ(5));
t187 = t150 * t164 + t157 * t168;
t186 = t150 * t168 - t157 * t164;
t183 = t162 * t165;
t182 = t162 * t169;
t180 = t164 * t165;
t179 = t164 * t170;
t178 = t165 * t168;
t175 = t168 * t170;
t172 = -t158 * t169 + t165 * t181;
t160 = -t163 * t177 + t173;
t159 = t163 * t176 + t174;
t156 = t163 * t165 + t166 * t182;
t155 = -t163 * t169 + t166 * t183;
t154 = t160 * t169 + t167 * t183;
t153 = t160 * t165 - t167 * t182;
t149 = t153 * t164 + t159 * t168;
t148 = -t153 * t168 + t159 * t164;
t1 = [-t187, -t159 * t180 + t160 * t168, t154 * t164, 0, -t148, 0; t149, -t157 * t180 + t158 * t168, -t172 * t164, 0, t186, 0; 0 (t165 * t179 + t166 * t168) * t162, t156 * t164, 0, t155 * t168 + t162 * t179, 0; t172, -t159 * t169, -t153, 0, 0, 0; t154, -t157 * t169, -t150, 0, 0, 0; 0, t170 * t182, -t155, 0, 0, 0; t186, t159 * t178 + t160 * t164, -t154 * t168, 0, t149, 0; t148, t157 * t178 + t158 * t164, t172 * t168, 0, t187, 0; 0 (t164 * t166 - t165 * t175) * t162, -t156 * t168, 0, t155 * t164 - t162 * t175, 0;];
JR_rot  = t1;
