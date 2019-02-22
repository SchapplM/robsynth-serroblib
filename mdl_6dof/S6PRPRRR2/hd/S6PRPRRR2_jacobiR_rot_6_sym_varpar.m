% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:38
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRR2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:38:25
% EndTime: 2019-02-22 09:38:25
% DurationCPUTime: 0.15s
% Computational Cost: add. (196->29), mult. (427->65), div. (0->0), fcn. (608->12), ass. (0->40)
t172 = sin(pkin(12));
t175 = cos(pkin(12));
t179 = sin(qJ(2));
t181 = cos(qJ(2));
t165 = t179 * t172 - t181 * t175;
t171 = qJ(5) + qJ(6);
t169 = sin(t171);
t180 = cos(qJ(4));
t191 = t169 * t180;
t170 = cos(t171);
t190 = t170 * t180;
t174 = sin(pkin(6));
t178 = sin(qJ(4));
t189 = t174 * t178;
t188 = t174 * t180;
t177 = cos(pkin(6));
t183 = t181 * t172 + t179 * t175;
t164 = t183 * t177;
t173 = sin(pkin(11));
t176 = cos(pkin(11));
t185 = t176 * t164 - t173 * t165;
t184 = -t173 * t164 - t176 * t165;
t182 = t165 * t177;
t163 = t183 * t174;
t162 = t165 * t174;
t160 = t163 * t180 + t177 * t178;
t159 = -t163 * t178 + t177 * t180;
t157 = t173 * t182 - t176 * t183;
t154 = -t173 * t183 - t176 * t182;
t152 = t173 * t189 + t180 * t184;
t151 = t173 * t188 - t178 * t184;
t150 = -t176 * t189 + t180 * t185;
t149 = -t176 * t188 - t178 * t185;
t148 = -t160 * t170 - t162 * t169;
t147 = -t160 * t169 + t162 * t170;
t146 = -t152 * t170 + t157 * t169;
t145 = -t152 * t169 - t157 * t170;
t144 = -t150 * t170 + t154 * t169;
t143 = -t150 * t169 - t154 * t170;
t1 = [0, t157 * t190 + t169 * t184, 0, t151 * t170, t145, t145; 0, t154 * t190 + t169 * t185, 0, t149 * t170, t143, t143; 0, -t162 * t190 + t163 * t169, 0, t159 * t170, t147, t147; 0, -t157 * t191 + t170 * t184, 0, -t151 * t169, t146, t146; 0, -t154 * t191 + t170 * t185, 0, -t149 * t169, t144, t144; 0, t162 * t191 + t163 * t170, 0, -t159 * t169, t148, t148; 0, t157 * t178, 0, t152, 0, 0; 0, t154 * t178, 0, t150, 0, 0; 0, -t162 * t178, 0, t160, 0, 0;];
JR_rot  = t1;
