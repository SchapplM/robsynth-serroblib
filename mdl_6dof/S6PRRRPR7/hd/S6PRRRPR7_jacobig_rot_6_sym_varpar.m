% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPR7_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:04
% EndTime: 2019-02-26 20:14:04
% DurationCPUTime: 0.10s
% Computational Cost: add. (34->21), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->29)
t170 = sin(pkin(12));
t172 = sin(pkin(6));
t190 = t170 * t172;
t171 = sin(pkin(7));
t189 = t171 * t172;
t175 = cos(pkin(6));
t188 = t171 * t175;
t173 = cos(pkin(12));
t187 = t173 * t172;
t174 = cos(pkin(7));
t181 = cos(qJ(2));
t186 = t174 * t181;
t178 = sin(qJ(2));
t185 = t175 * t178;
t184 = t175 * t181;
t166 = -t170 * t178 + t173 * t184;
t183 = -t166 * t174 + t171 * t187;
t168 = -t170 * t184 - t173 * t178;
t182 = t168 * t174 + t170 * t189;
t180 = cos(qJ(3));
t179 = cos(qJ(4));
t177 = sin(qJ(3));
t176 = sin(qJ(4));
t169 = -t170 * t185 + t173 * t181;
t167 = t170 * t181 + t173 * t185;
t165 = t175 * t174 - t181 * t189;
t164 = -t168 * t171 + t174 * t190;
t163 = -t166 * t171 - t174 * t187;
t1 = [0, t190, t164, t169 * t177 - t182 * t180, 0 (t169 * t180 + t182 * t177) * t176 - t164 * t179; 0, -t187, t163, t167 * t177 + t183 * t180, 0 (t167 * t180 - t183 * t177) * t176 - t163 * t179; 0, t175, t165, -t180 * t188 + (t177 * t178 - t180 * t186) * t172, 0 (t177 * t188 + (t177 * t186 + t178 * t180) * t172) * t176 - t165 * t179;];
Jg_rot  = t1;
