% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRP8_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:43:49
% EndTime: 2019-02-26 22:43:49
% DurationCPUTime: 0.03s
% Computational Cost: add. (18->11), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
t173 = sin(pkin(6));
t177 = cos(qJ(2));
t185 = t173 * t177;
t176 = sin(qJ(1));
t184 = t176 * t173;
t175 = sin(qJ(2));
t183 = t176 * t175;
t182 = t176 * t177;
t178 = cos(qJ(1));
t181 = t178 * t173;
t180 = t178 * t175;
t179 = t178 * t177;
t174 = cos(pkin(6));
t172 = qJ(3) + qJ(4);
t171 = cos(t172);
t170 = sin(t172);
t169 = t174 * t182 + t180;
t168 = -t174 * t179 + t183;
t1 = [0, t184, t169, t169 (-t174 * t183 + t179) * t170 - t171 * t184, 0; 0, -t181, t168, t168 (t174 * t180 + t182) * t170 + t171 * t181, 0; 1, t174, -t185, -t185, t173 * t175 * t170 - t174 * t171, 0;];
Jg_rot  = t1;
