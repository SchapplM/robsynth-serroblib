% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPR8_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:37
% EndTime: 2019-02-26 20:14:37
% DurationCPUTime: 0.06s
% Computational Cost: add. (34->21), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->29)
t164 = sin(pkin(12));
t166 = sin(pkin(6));
t184 = t164 * t166;
t165 = sin(pkin(7));
t183 = t165 * t166;
t169 = cos(pkin(6));
t182 = t165 * t169;
t167 = cos(pkin(12));
t181 = t167 * t166;
t168 = cos(pkin(7));
t175 = cos(qJ(2));
t180 = t168 * t175;
t172 = sin(qJ(2));
t179 = t169 * t172;
t178 = t169 * t175;
t160 = -t164 * t172 + t167 * t178;
t177 = -t160 * t168 + t165 * t181;
t162 = -t164 * t178 - t167 * t172;
t176 = t162 * t168 + t164 * t183;
t174 = cos(qJ(3));
t173 = cos(qJ(4));
t171 = sin(qJ(3));
t170 = sin(qJ(4));
t163 = -t164 * t179 + t167 * t175;
t161 = t164 * t175 + t167 * t179;
t159 = t169 * t168 - t175 * t183;
t158 = -t162 * t165 + t168 * t184;
t157 = -t160 * t165 - t168 * t181;
t1 = [0, t184, t158, t163 * t171 - t176 * t174, 0 (t163 * t174 + t176 * t171) * t173 + t158 * t170; 0, -t181, t157, t161 * t171 + t177 * t174, 0 (t161 * t174 - t177 * t171) * t173 + t157 * t170; 0, t169, t159, -t174 * t182 + (t171 * t172 - t174 * t180) * t166, 0 (t171 * t182 + (t171 * t180 + t172 * t174) * t166) * t173 + t159 * t170;];
Jg_rot  = t1;
