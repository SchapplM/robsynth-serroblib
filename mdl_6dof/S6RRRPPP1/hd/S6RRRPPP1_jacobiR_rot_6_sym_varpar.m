% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPP1_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:02:48
% EndTime: 2019-02-26 22:02:48
% DurationCPUTime: 0.09s
% Computational Cost: add. (53->27), mult. (183->71), div. (0->0), fcn. (258->10), ass. (0->35)
t159 = sin(pkin(10));
t162 = cos(pkin(6));
t188 = t159 * t162;
t160 = sin(pkin(6));
t164 = sin(qJ(2));
t187 = t160 * t164;
t161 = cos(pkin(10));
t186 = t161 * t162;
t166 = cos(qJ(3));
t185 = t162 * t166;
t163 = sin(qJ(3));
t167 = cos(qJ(2));
t184 = t163 * t167;
t183 = t164 * t162;
t165 = sin(qJ(1));
t182 = t164 * t165;
t181 = t164 * t166;
t168 = cos(qJ(1));
t180 = t164 * t168;
t179 = t165 * t163;
t178 = t166 * t167;
t177 = t168 * t163;
t176 = t168 * t166;
t155 = t167 * t179 + t176;
t175 = -t155 * t162 + t160 * t182;
t157 = t165 * t166 - t167 * t177;
t174 = t157 * t162 + t160 * t180;
t173 = t162 * t184 - t187;
t172 = t160 * t167 + t163 * t183;
t171 = t162 * t167 - t163 * t187;
t170 = -t159 * t181 - t172 * t161;
t169 = t172 * t159 - t161 * t181;
t158 = t167 * t176 + t179;
t156 = -t165 * t178 + t177;
t1 = [-t155 * t160 - t162 * t182, t171 * t168, t158 * t160, 0, 0, 0; -t157 * t160 + t162 * t180, t171 * t165, -t156 * t160, 0, 0, 0; 0, t160 * t184 + t183, t160 * t181, 0, 0, 0; t156 * t159 + t175 * t161, t170 * t168, t157 * t159 + t158 * t186, 0, 0, 0; t158 * t159 - t174 * t161, t170 * t165, -t155 * t159 - t156 * t186, 0, 0, 0; 0, t159 * t178 + t173 * t161 (-t159 * t163 + t161 * t185) * t164, 0, 0, 0; t156 * t161 - t175 * t159, t169 * t168, t157 * t161 - t158 * t188, 0, 0, 0; t158 * t161 + t174 * t159, t169 * t165, -t155 * t161 + t156 * t188, 0, 0, 0; 0, -t173 * t159 + t161 * t178 (-t159 * t185 - t161 * t163) * t164, 0, 0, 0;];
JR_rot  = t1;
