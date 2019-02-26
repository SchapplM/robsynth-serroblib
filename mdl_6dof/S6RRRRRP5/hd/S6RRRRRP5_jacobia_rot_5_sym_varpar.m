% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRP5_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:41:48
% EndTime: 2019-02-26 22:41:48
% DurationCPUTime: 0.11s
% Computational Cost: add. (281->21), mult. (278->55), div. (62->9), fcn. (402->9), ass. (0->35)
t63 = cos(qJ(2));
t61 = sin(qJ(2));
t62 = sin(qJ(1));
t69 = t62 * t61;
t53 = atan2(-t69, -t63);
t51 = sin(t53);
t52 = cos(t53);
t45 = -t51 * t69 - t52 * t63;
t44 = 0.1e1 / t45 ^ 2;
t64 = cos(qJ(1));
t74 = t44 * t64 ^ 2;
t57 = qJ(3) + qJ(4) + qJ(5);
t55 = sin(t57);
t56 = cos(t57);
t66 = t64 * t56;
t50 = t62 * t55 + t63 * t66;
t48 = 0.1e1 / t50 ^ 2;
t67 = t64 * t55;
t49 = -t62 * t56 + t63 * t67;
t73 = t48 * t49;
t58 = t61 ^ 2;
t72 = t58 / t63 ^ 2;
t71 = t61 * t64;
t54 = 0.1e1 / (t62 ^ 2 * t72 + 0.1e1);
t70 = t62 * t54;
t68 = t62 * t63;
t65 = t49 ^ 2 * t48 + 0.1e1;
t59 = 0.1e1 / t63;
t47 = 0.1e1 / t50;
t46 = (0.1e1 + t72) * t70;
t43 = 0.1e1 / t45;
t42 = 0.1e1 / (t58 * t74 + 0.1e1);
t41 = 0.1e1 / t65;
t40 = t65 * t41;
t1 = [t59 * t54 * t71, t46, 0, 0, 0, 0; (-t43 * t69 - (-t52 * t58 * t59 * t70 + (t54 - 0.1e1) * t61 * t51) * t61 * t74) * t42 (t63 * t43 - (-t51 * t68 + t52 * t61 + (t51 * t63 - t52 * t69) * t46) * t61 * t44) * t64 * t42, 0, 0, 0, 0; ((-t55 * t68 - t66) * t47 - (-t56 * t68 + t67) * t73) * t41 (-t47 * t55 + t56 * t73) * t41 * t71, t40, t40, t40, 0;];
Ja_rot  = t1;
