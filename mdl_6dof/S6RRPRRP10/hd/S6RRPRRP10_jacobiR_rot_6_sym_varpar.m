% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP10_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:51:04
% EndTime: 2019-02-26 21:51:04
% DurationCPUTime: 0.10s
% Computational Cost: add. (122->30), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
t174 = cos(pkin(6));
t176 = sin(qJ(2));
t180 = cos(qJ(1));
t182 = t180 * t176;
t177 = sin(qJ(1));
t179 = cos(qJ(2));
t184 = t177 * t179;
t165 = t174 * t182 + t184;
t172 = pkin(11) + qJ(4);
t170 = sin(t172);
t171 = cos(t172);
t173 = sin(pkin(6));
t187 = t173 * t180;
t159 = -t165 * t171 + t170 * t187;
t181 = t180 * t179;
t185 = t177 * t176;
t164 = -t174 * t181 + t185;
t175 = sin(qJ(5));
t178 = cos(qJ(5));
t195 = t159 * t175 + t164 * t178;
t194 = t159 * t178 - t164 * t175;
t191 = t171 * t175;
t190 = t171 * t178;
t189 = t173 * t176;
t188 = t173 * t177;
t186 = t175 * t179;
t183 = t178 * t179;
t157 = -t165 * t170 - t171 * t187;
t167 = -t174 * t185 + t181;
t166 = t174 * t184 + t182;
t163 = t174 * t170 + t171 * t189;
t162 = -t170 * t189 + t174 * t171;
t161 = t167 * t171 + t170 * t188;
t160 = t167 * t170 - t171 * t188;
t156 = t161 * t178 + t166 * t175;
t155 = t161 * t175 - t166 * t178;
t1 = [t194, -t166 * t190 + t167 * t175, 0, -t160 * t178, -t155, 0; t156, -t164 * t190 + t165 * t175, 0, t157 * t178, t195, 0; 0 (t171 * t183 + t175 * t176) * t173, 0, t162 * t178, -t163 * t175 - t173 * t183, 0; t157, -t166 * t170, 0, t161, 0, 0; t160, -t164 * t170, 0, -t159, 0, 0; 0, t173 * t179 * t170, 0, t163, 0, 0; t195, -t166 * t191 - t167 * t178, 0, -t160 * t175, t156, 0; t155, -t164 * t191 - t165 * t178, 0, t157 * t175, -t194, 0; 0 (t171 * t186 - t176 * t178) * t173, 0, t162 * t175, t163 * t178 - t173 * t186, 0;];
JR_rot  = t1;
