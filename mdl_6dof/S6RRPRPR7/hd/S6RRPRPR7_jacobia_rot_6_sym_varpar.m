% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRPRPR7_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:12
% EndTime: 2019-02-26 21:41:13
% DurationCPUTime: 0.15s
% Computational Cost: add. (737->24), mult. (876->55), div. (85->9), fcn. (1282->11), ass. (0->38)
t87 = qJ(4) + pkin(10);
t73 = sin(t87);
t77 = cos(qJ(2));
t84 = cos(t87);
t95 = sin(qJ(2));
t67 = -t77 * t73 + t95 * t84;
t75 = sin(qJ(1));
t59 = t67 * t75;
t66 = t95 * t73 + t77 * t84;
t54 = atan2(t59, t66);
t51 = sin(t54);
t52 = cos(t54);
t49 = t51 * t59 + t52 * t66;
t48 = 0.1e1 / t49 ^ 2;
t96 = cos(qJ(1));
t63 = t67 * t96;
t89 = t63 ^ 2 * t48;
t46 = 0.1e1 / (0.1e1 + t89);
t47 = 0.1e1 / t49;
t61 = t66 * t75;
t62 = t66 * t96;
t92 = t52 * t59;
t65 = 0.1e1 / t66 ^ 2;
t53 = 0.1e1 / (t59 ^ 2 * t65 + 0.1e1);
t64 = 0.1e1 / t66;
t98 = (-t59 * t65 * t67 - t61 * t64) * t53;
t101 = ((-t51 * t61 + t52 * t67 + (-t51 * t66 + t92) * t98) * t48 * t63 + t62 * t47) * t46;
t74 = sin(qJ(6));
t76 = cos(qJ(6));
t58 = t62 * t76 - t75 * t74;
t56 = 0.1e1 / t58 ^ 2;
t57 = t62 * t74 + t75 * t76;
t86 = t57 ^ 2 * t56 + 0.1e1;
t50 = 0.1e1 / t86;
t55 = 0.1e1 / t58;
t91 = t56 * t57;
t97 = (-t74 * t55 + t76 * t91) * t50 * t63;
t1 = [t63 * t64 * t53, -t98, 0, t98, 0, 0; (t59 * t47 - (-t51 + (-t64 * t92 + t51) * t53) * t89) * t46, -t101, 0, t101, 0, 0; ((-t61 * t74 + t96 * t76) * t55 - (-t61 * t76 - t96 * t74) * t91) * t50, t97, 0, -t97, 0, t86 * t50;];
Ja_rot  = t1;
