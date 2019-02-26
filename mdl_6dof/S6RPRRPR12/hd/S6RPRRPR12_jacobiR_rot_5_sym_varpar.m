% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRRPR12_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:07:16
% EndTime: 2019-02-26 21:07:16
% DurationCPUTime: 0.12s
% Computational Cost: add. (105->30), mult. (307->59), div. (0->0), fcn. (430->12), ass. (0->40)
t180 = cos(pkin(6));
t175 = sin(pkin(12));
t186 = cos(qJ(1));
t190 = t186 * t175;
t178 = cos(pkin(12));
t183 = sin(qJ(1));
t191 = t183 * t178;
t170 = t180 * t190 + t191;
t182 = sin(qJ(3));
t185 = cos(qJ(3));
t189 = t186 * t178;
t192 = t183 * t175;
t169 = -t180 * t189 + t192;
t176 = sin(pkin(7));
t179 = cos(pkin(7));
t177 = sin(pkin(6));
t194 = t177 * t186;
t187 = t169 * t179 + t176 * t194;
t159 = -t170 * t185 + t187 * t182;
t164 = -t169 * t176 + t179 * t194;
t181 = sin(qJ(4));
t184 = cos(qJ(4));
t201 = t159 * t181 - t164 * t184;
t200 = -t159 * t184 - t164 * t181;
t196 = t176 * t180;
t195 = t177 * t183;
t193 = t179 * t185;
t188 = t176 * t195;
t157 = -t170 * t182 - t187 * t185;
t172 = -t180 * t192 + t189;
t171 = -t180 * t191 - t190;
t168 = -t177 * t178 * t176 + t180 * t179;
t166 = -t171 * t176 + t179 * t195;
t163 = t182 * t196 + (t178 * t179 * t182 + t175 * t185) * t177;
t162 = t185 * t196 + (-t175 * t182 + t178 * t193) * t177;
t161 = t172 * t185 + (t171 * t179 + t188) * t182;
t160 = -t171 * t193 + t172 * t182 - t185 * t188;
t156 = t161 * t184 + t166 * t181;
t155 = t161 * t181 - t166 * t184;
t1 = [t157, 0, t161, 0, 0, 0; t160, 0, -t159, 0, 0, 0; 0, 0, t163, 0, 0, 0; t200, 0, t160 * t184, t155, 0, 0; -t156, 0, -t157 * t184, -t201, 0, 0; 0, 0, -t162 * t184, t163 * t181 - t168 * t184, 0, 0; t201, 0, -t160 * t181, t156, 0, 0; t155, 0, t157 * t181, t200, 0, 0; 0, 0, t162 * t181, t163 * t184 + t168 * t181, 0, 0;];
JR_rot  = t1;
