% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRR10_jacobig_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobig_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobig_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:57
% EndTime: 2019-02-26 22:52:57
% DurationCPUTime: 0.05s
% Computational Cost: add. (24->17), mult. (69->37), div. (0->0), fcn. (101->12), ass. (0->24)
t175 = sin(pkin(6));
t181 = sin(qJ(1));
t190 = t181 * t175;
t180 = sin(qJ(2));
t189 = t181 * t180;
t183 = cos(qJ(2));
t188 = t181 * t183;
t184 = cos(qJ(1));
t187 = t184 * t175;
t186 = t184 * t180;
t185 = t184 * t183;
t182 = cos(qJ(3));
t179 = sin(qJ(3));
t178 = cos(pkin(6));
t177 = cos(pkin(7));
t176 = cos(pkin(8));
t174 = sin(pkin(7));
t173 = sin(pkin(8));
t172 = -t178 * t188 - t186;
t171 = t178 * t185 - t189;
t170 = -t175 * t183 * t174 + t178 * t177;
t169 = -t172 * t174 + t177 * t190;
t168 = -t171 * t174 - t177 * t187;
t1 = [0, t190, t169 -(-(-t178 * t189 + t185) * t179 + (t172 * t177 + t174 * t190) * t182) * t173 + t169 * t176, 0, 0; 0, -t187, t168 -(-(t178 * t186 + t188) * t179 + (t171 * t177 - t174 * t187) * t182) * t173 + t168 * t176, 0, 0; 1, t178, t170 -(t178 * t174 * t182 + (t177 * t182 * t183 - t179 * t180) * t175) * t173 + t170 * t176, 0, 0;];
Jg_rot  = t1;
