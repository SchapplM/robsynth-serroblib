% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR7_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobig_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:19
% EndTime: 2019-02-26 22:19:19
% DurationCPUTime: 0.03s
% Computational Cost: add. (24->11), mult. (31->18), div. (0->0), fcn. (52->8), ass. (0->19)
t133 = sin(pkin(6));
t137 = cos(qJ(2));
t145 = t133 * t137;
t136 = sin(qJ(1));
t144 = t136 * t133;
t135 = sin(qJ(2));
t143 = t136 * t135;
t142 = t136 * t137;
t138 = cos(qJ(1));
t141 = t138 * t133;
t140 = t138 * t135;
t139 = t138 * t137;
t134 = cos(pkin(6));
t132 = qJ(3) + pkin(12) + qJ(5);
t131 = cos(t132);
t130 = sin(t132);
t129 = t134 * t142 + t140;
t128 = -t134 * t139 + t143;
t1 = [0, t144, t129, 0, t129 (-t134 * t143 + t139) * t130 - t131 * t144; 0, -t141, t128, 0, t128 (t134 * t140 + t142) * t130 + t131 * t141; 1, t134, -t145, 0, -t145, t133 * t135 * t130 - t134 * t131;];
Jg_rot  = t1;
