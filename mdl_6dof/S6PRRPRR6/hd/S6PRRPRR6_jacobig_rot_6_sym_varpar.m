% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR6_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:06
% EndTime: 2019-02-26 20:07:06
% DurationCPUTime: 0.08s
% Computational Cost: add. (40->22), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->30)
t171 = sin(pkin(12));
t173 = sin(pkin(6));
t189 = t171 * t173;
t172 = sin(pkin(7));
t188 = t172 * t173;
t176 = cos(pkin(6));
t187 = t172 * t176;
t174 = cos(pkin(12));
t186 = t174 * t173;
t175 = cos(pkin(7));
t180 = cos(qJ(2));
t185 = t175 * t180;
t178 = sin(qJ(2));
t184 = t176 * t178;
t183 = t176 * t180;
t164 = -t171 * t178 + t174 * t183;
t182 = -t164 * t175 + t172 * t186;
t166 = -t171 * t183 - t174 * t178;
t181 = t166 * t175 + t171 * t188;
t179 = cos(qJ(3));
t177 = sin(qJ(3));
t170 = pkin(13) + qJ(5);
t169 = cos(t170);
t168 = sin(t170);
t167 = -t171 * t184 + t174 * t180;
t165 = t171 * t180 + t174 * t184;
t163 = t176 * t175 - t180 * t188;
t162 = -t166 * t172 + t175 * t189;
t161 = -t164 * t172 - t175 * t186;
t1 = [0, t189, t162, 0, t167 * t177 - t181 * t179 (t167 * t179 + t181 * t177) * t168 - t162 * t169; 0, -t186, t161, 0, t165 * t177 + t182 * t179 (t165 * t179 - t182 * t177) * t168 - t161 * t169; 0, t176, t163, 0, -t179 * t187 + (t177 * t178 - t179 * t185) * t173 (t177 * t187 + (t177 * t185 + t178 * t179) * t173) * t168 - t163 * t169;];
Jg_rot  = t1;
