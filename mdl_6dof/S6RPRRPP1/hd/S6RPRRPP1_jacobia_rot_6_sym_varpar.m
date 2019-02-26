% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RPRRPP1_jacobia_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_jacobia_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_jacobia_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:56:16
% EndTime: 2019-02-26 20:56:17
% DurationCPUTime: 0.14s
% Computational Cost: add. (665->28), mult. (509->70), div. (100->11), fcn. (770->9), ass. (0->40)
t65 = sin(qJ(3));
t77 = t65 ^ 2;
t60 = qJ(4) + pkin(10);
t56 = sin(t60);
t58 = cos(t60);
t61 = qJ(1) + pkin(9);
t59 = cos(t61);
t57 = sin(t61);
t66 = cos(qJ(3));
t71 = t57 * t66;
t46 = t56 * t71 + t59 * t58;
t68 = t65 * t56;
t43 = atan2(-t46, t68);
t39 = sin(t43);
t40 = cos(t43);
t38 = -t39 * t46 + t40 * t68;
t37 = 0.1e1 / t38 ^ 2;
t69 = t59 * t66;
t49 = t56 * t69 - t57 * t58;
t76 = t37 * t49;
t74 = t40 * t46;
t73 = t49 ^ 2 * t37;
t53 = 0.1e1 / t56;
t63 = 0.1e1 / t65;
t72 = t53 * t63;
t70 = t59 * t65;
t50 = t57 * t56 + t58 * t69;
t45 = 0.1e1 / t50 ^ 2;
t67 = t59 ^ 2 * t77 * t45;
t64 = 0.1e1 / t77;
t54 = 0.1e1 / t56 ^ 2;
t48 = -t59 * t56 + t58 * t71;
t44 = 0.1e1 / t50;
t42 = 0.1e1 / (t46 ^ 2 * t64 * t54 + 0.1e1);
t41 = 0.1e1 / (0.1e1 + t67);
t36 = 0.1e1 / t38;
t35 = (t46 * t53 * t64 * t66 + t57) * t42;
t34 = 0.1e1 / (0.1e1 + t73);
t33 = (t46 * t54 * t58 - t48 * t53) * t63 * t42;
t1 = [-t49 * t42 * t72, 0, t35, t33, 0, 0; (-t46 * t36 - (-t39 + (t72 * t74 + t39) * t42) * t73) * t34, 0 (t35 * t74 * t76 + (-t36 * t70 - (t40 * t66 + (-t35 + t57) * t65 * t39) * t76) * t56) * t34 (t50 * t36 - (t40 * t65 * t58 - t39 * t48 + (-t39 * t68 - t74) * t33) * t76) * t34, 0, 0; (-t45 * t48 * t59 + t44 * t57) * t65 * t41, 0 (-t44 * t69 - t58 * t67) * t41, -t49 * t45 * t41 * t70, 0, 0;];
Ja_rot  = t1;
