% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:25
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR5_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:25:09
% EndTime: 2019-02-22 11:25:09
% DurationCPUTime: 0.12s
% Computational Cost: add. (192->32), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->38)
t176 = sin(qJ(4));
t179 = cos(qJ(4));
t175 = cos(pkin(6));
t172 = sin(pkin(11));
t174 = cos(pkin(11));
t177 = sin(qJ(2));
t180 = cos(qJ(2));
t183 = t180 * t172 + t177 * t174;
t163 = t183 * t175;
t164 = t177 * t172 - t180 * t174;
t178 = sin(qJ(1));
t181 = cos(qJ(1));
t185 = t181 * t163 - t178 * t164;
t173 = sin(pkin(6));
t188 = t173 * t181;
t149 = t176 * t188 - t179 * t185;
t182 = t164 * t175;
t153 = -t178 * t183 - t181 * t182;
t171 = pkin(12) + qJ(6);
t169 = sin(t171);
t170 = cos(t171);
t195 = t149 * t169 - t153 * t170;
t194 = t149 * t170 + t153 * t169;
t191 = t169 * t179;
t190 = t170 * t179;
t189 = t173 * t178;
t184 = -t178 * t163 - t181 * t164;
t147 = -t176 * t185 - t179 * t188;
t162 = t183 * t173;
t161 = t164 * t173;
t159 = t162 * t179 + t175 * t176;
t158 = -t162 * t176 + t175 * t179;
t156 = t178 * t182 - t181 * t183;
t151 = t176 * t189 + t179 * t184;
t150 = t176 * t184 - t179 * t189;
t146 = t151 * t170 - t156 * t169;
t145 = -t151 * t169 - t156 * t170;
t1 = [t194, t156 * t190 + t169 * t184, 0, -t150 * t170, 0, t145; t146, t153 * t190 + t169 * t185, 0, t147 * t170, 0, t195; 0, -t161 * t190 + t162 * t169, 0, t158 * t170, 0, -t159 * t169 + t161 * t170; -t195, -t156 * t191 + t170 * t184, 0, t150 * t169, 0, -t146; t145, -t153 * t191 + t170 * t185, 0, -t147 * t169, 0, t194; 0, t161 * t191 + t162 * t170, 0, -t158 * t169, 0, -t159 * t170 - t161 * t169; t147, t156 * t176, 0, t151, 0, 0; t150, t153 * t176, 0, -t149, 0, 0; 0, -t161 * t176, 0, t159, 0, 0;];
JR_rot  = t1;
