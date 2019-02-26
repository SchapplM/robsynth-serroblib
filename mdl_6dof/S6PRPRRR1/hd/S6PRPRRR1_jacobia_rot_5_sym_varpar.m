% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPRRR1_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobia_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:43
% EndTime: 2019-02-26 19:53:43
% DurationCPUTime: 0.09s
% Computational Cost: add. (297->19), mult. (709->46), div. (40->9), fcn. (992->13), ass. (0->36)
t74 = sin(pkin(11));
t75 = sin(pkin(6));
t84 = t74 * t75;
t78 = cos(pkin(6));
t73 = sin(pkin(12));
t76 = cos(pkin(12));
t79 = sin(qJ(2));
t80 = cos(qJ(2));
t82 = t80 * t73 + t79 * t76;
t67 = t82 * t78;
t68 = t79 * t73 - t80 * t76;
t77 = cos(pkin(11));
t62 = -t74 * t67 - t77 * t68;
t72 = qJ(4) + qJ(5);
t70 = sin(t72);
t71 = cos(t72);
t57 = t62 * t71 + t70 * t84;
t55 = 0.1e1 / t57 ^ 2;
t56 = t62 * t70 - t71 * t84;
t83 = t56 ^ 2 * t55 + 0.1e1;
t81 = t68 * t78;
t66 = t82 * t75;
t65 = t68 * t75;
t64 = 0.1e1 / t65 ^ 2;
t60 = t74 * t81 - t77 * t82;
t59 = -t74 * t82 - t77 * t81;
t58 = -t77 * t67 + t74 * t68;
t54 = atan2(t59, t65);
t52 = cos(t54);
t51 = sin(t54);
t50 = 0.1e1 / t83;
t49 = t51 * t59 + t52 * t65;
t48 = 0.1e1 / t49 ^ 2;
t46 = (t58 / t65 - t66 * t59 * t64) / (t59 ^ 2 * t64 + 0.1e1);
t45 = t83 * t50;
t1 = [0, t46, 0, 0, 0, 0; 0 (t62 / t49 + (t51 * t58 + t52 * t66 + (-t51 * t65 + t52 * t59) * t46) * t60 * t48) / (t60 ^ 2 * t48 + 0.1e1) 0, 0, 0, 0; 0 (t70 / t57 - t71 * t56 * t55) * t60 * t50, 0, t45, t45, 0;];
Ja_rot  = t1;
