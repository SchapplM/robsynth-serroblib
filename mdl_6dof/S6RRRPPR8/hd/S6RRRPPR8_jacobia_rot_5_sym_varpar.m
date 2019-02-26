% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRPPR8_jacobia_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobia_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobia_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:39
% EndTime: 2019-02-26 22:07:39
% DurationCPUTime: 0.09s
% Computational Cost: add. (151->23), mult. (457->62), div. (73->11), fcn. (682->11), ass. (0->42)
t61 = cos(pkin(6));
t66 = cos(qJ(2));
t67 = cos(qJ(1));
t68 = t67 * t66;
t63 = sin(qJ(2));
t64 = sin(qJ(1));
t71 = t64 * t63;
t55 = -t61 * t71 + t68;
t62 = sin(qJ(3));
t65 = cos(qJ(3));
t60 = sin(pkin(6));
t74 = t60 * t64;
t45 = t55 * t62 - t65 * t74;
t43 = 0.1e1 / t45 ^ 2;
t46 = t55 * t65 + t62 * t74;
t79 = t43 * t46;
t78 = t46 ^ 2 * t43;
t51 = -t61 * t68 + t71;
t73 = t60 * t66;
t50 = atan2(t51, t73);
t48 = cos(t50);
t77 = t48 * t51;
t47 = sin(t50);
t40 = t47 * t51 + t48 * t73;
t39 = 0.1e1 / t40 ^ 2;
t69 = t67 * t63;
t70 = t64 * t66;
t53 = t61 * t70 + t69;
t76 = t53 ^ 2 * t39;
t57 = 0.1e1 / t60;
t58 = 0.1e1 / t66;
t75 = t57 * t58;
t72 = t60 * t67;
t59 = 0.1e1 / t66 ^ 2;
t52 = t61 * t69 + t70;
t49 = 0.1e1 / (0.1e1 + t51 ^ 2 / t60 ^ 2 * t59);
t42 = 0.1e1 / t45;
t41 = 0.1e1 / (0.1e1 + t78);
t38 = 0.1e1 / t40;
t37 = 0.1e1 / (0.1e1 + t76);
t36 = (t51 * t59 * t63 + t52 * t58) * t57 * t49;
t1 = [t53 * t49 * t75, t36, 0, 0, 0, 0; (t51 * t38 + (t47 + (t75 * t77 - t47) * t49) * t76) * t37 (-t55 * t38 + (-t48 * t60 * t63 + t47 * t52 + (-t47 * t73 + t77) * t36) * t53 * t39) * t37, 0, 0, 0, 0; ((-t52 * t65 + t62 * t72) * t42 - (-t52 * t62 - t65 * t72) * t79) * t41 (-t65 * t42 + t62 * t79) * t53 * t41 (-t42 * t45 - t78) * t41, 0, 0, 0;];
Ja_rot  = t1;
