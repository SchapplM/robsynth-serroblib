% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6PRPPRR3_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobia_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:50
% EndTime: 2019-02-26 19:45:50
% DurationCPUTime: 0.11s
% Computational Cost: add. (222->23), mult. (626->54), div. (35->9), fcn. (880->13), ass. (0->36)
t60 = sin(pkin(10));
t61 = sin(pkin(6));
t73 = t60 * t61;
t64 = cos(pkin(6));
t66 = sin(qJ(2));
t72 = t64 * t66;
t68 = cos(qJ(2));
t71 = t64 * t68;
t63 = cos(pkin(10));
t57 = t60 * t71 + t63 * t66;
t58 = -t60 * t72 + t63 * t68;
t59 = sin(pkin(11));
t62 = cos(pkin(11));
t50 = t57 * t59 + t58 * t62;
t65 = sin(qJ(5));
t67 = cos(qJ(5));
t44 = t50 * t67 - t65 * t73;
t42 = 0.1e1 / t44 ^ 2;
t43 = t50 * t65 + t67 * t73;
t70 = t43 ^ 2 * t42 + 0.1e1;
t49 = -t57 * t62 + t58 * t59;
t69 = -t60 * t66 + t63 * t71;
t56 = t60 * t68 + t63 * t72;
t55 = (t59 * t68 - t62 * t66) * t61;
t54 = (t59 * t66 + t62 * t68) * t61;
t53 = 0.1e1 / t54 ^ 2;
t46 = t56 * t59 + t69 * t62;
t45 = -t56 * t62 + t69 * t59;
t41 = atan2(-t46, t54);
t39 = cos(t41);
t38 = sin(t41);
t37 = 0.1e1 / t70;
t36 = -t38 * t46 + t39 * t54;
t35 = 0.1e1 / t36 ^ 2;
t33 = (-t45 / t54 + t55 * t46 * t53) / (t46 ^ 2 * t53 + 0.1e1);
t1 = [0, t33, 0, 0, 0, 0; 0 (-t50 / t36 - (-t38 * t45 + t39 * t55 + (-t38 * t54 - t39 * t46) * t33) * t49 * t35) / (t49 ^ 2 * t35 + 0.1e1) 0, 0, 0, 0; 0 (t65 / t44 - t67 * t43 * t42) * t49 * t37, 0, 0, t70 * t37, 0;];
Ja_rot  = t1;
