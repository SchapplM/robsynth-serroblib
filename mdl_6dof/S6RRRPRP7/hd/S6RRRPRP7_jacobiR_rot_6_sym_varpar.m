% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRP7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:39
% EndTime: 2019-02-26 22:12:39
% DurationCPUTime: 0.11s
% Computational Cost: add. (122->30), mult. (219->63), div. (0->0), fcn. (320->10), ass. (0->37)
t175 = cos(pkin(6));
t177 = sin(qJ(2));
t181 = cos(qJ(1));
t183 = t181 * t177;
t178 = sin(qJ(1));
t180 = cos(qJ(2));
t185 = t178 * t180;
t166 = t175 * t183 + t185;
t173 = qJ(3) + pkin(11);
t171 = sin(t173);
t172 = cos(t173);
t174 = sin(pkin(6));
t188 = t174 * t181;
t160 = -t166 * t172 + t171 * t188;
t182 = t181 * t180;
t186 = t178 * t177;
t165 = -t175 * t182 + t186;
t176 = sin(qJ(5));
t179 = cos(qJ(5));
t196 = t160 * t176 + t165 * t179;
t195 = t160 * t179 - t165 * t176;
t192 = t172 * t176;
t191 = t172 * t179;
t190 = t174 * t177;
t189 = t174 * t178;
t187 = t176 * t180;
t184 = t179 * t180;
t158 = -t166 * t171 - t172 * t188;
t168 = -t175 * t186 + t182;
t167 = t175 * t185 + t183;
t164 = t175 * t171 + t172 * t190;
t163 = -t171 * t190 + t175 * t172;
t162 = t168 * t172 + t171 * t189;
t161 = t168 * t171 - t172 * t189;
t157 = t162 * t179 + t167 * t176;
t156 = t162 * t176 - t167 * t179;
t1 = [t195, -t167 * t191 + t168 * t176, -t161 * t179, 0, -t156, 0; t157, -t165 * t191 + t166 * t176, t158 * t179, 0, t196, 0; 0 (t172 * t184 + t176 * t177) * t174, t163 * t179, 0, -t164 * t176 - t174 * t184, 0; t158, -t167 * t171, t162, 0, 0, 0; t161, -t165 * t171, -t160, 0, 0, 0; 0, t174 * t180 * t171, t164, 0, 0, 0; t196, -t167 * t192 - t168 * t179, -t161 * t176, 0, t157, 0; t156, -t165 * t192 - t166 * t179, t158 * t176, 0, -t195, 0; 0 (t172 * t187 - t177 * t179) * t174, t163 * t176, 0, t164 * t179 - t174 * t187, 0;];
JR_rot  = t1;
