% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRRP5_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:17:29
% EndTime: 2019-02-26 20:17:29
% DurationCPUTime: 0.06s
% Computational Cost: add. (34->21), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->29)
t171 = sin(pkin(12));
t173 = sin(pkin(6));
t191 = t171 * t173;
t172 = sin(pkin(7));
t190 = t172 * t173;
t176 = cos(pkin(6));
t189 = t172 * t176;
t174 = cos(pkin(12));
t188 = t174 * t173;
t175 = cos(pkin(7));
t182 = cos(qJ(2));
t187 = t175 * t182;
t179 = sin(qJ(2));
t186 = t176 * t179;
t185 = t176 * t182;
t167 = -t171 * t179 + t174 * t185;
t184 = -t167 * t175 + t172 * t188;
t169 = -t171 * t185 - t174 * t179;
t183 = t169 * t175 + t171 * t190;
t181 = cos(qJ(3));
t180 = cos(qJ(4));
t178 = sin(qJ(3));
t177 = sin(qJ(4));
t170 = -t171 * t186 + t174 * t182;
t168 = t171 * t182 + t174 * t186;
t166 = t176 * t175 - t182 * t190;
t165 = -t169 * t172 + t175 * t191;
t164 = -t167 * t172 - t175 * t188;
t1 = [0, t191, t165, t170 * t178 - t183 * t181 (t170 * t181 + t183 * t178) * t177 - t165 * t180, 0; 0, -t188, t164, t168 * t178 + t184 * t181 (t168 * t181 - t184 * t178) * t177 - t164 * t180, 0; 0, t176, t166, -t181 * t189 + (t178 * t179 - t181 * t187) * t173 (t178 * t189 + (t178 * t187 + t179 * t181) * t173) * t177 - t166 * t180, 0;];
Jg_rot  = t1;
